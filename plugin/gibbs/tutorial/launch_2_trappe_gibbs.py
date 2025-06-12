"""
Example Gibbs ensemble Monte Carlo simulation of TraPPE alkanes using FEASST.
Run multiple simulations, each on a single processor, using the argument "dictionary_input."
In particular, list the temperatures to run for each particle, as well as the expected densities and number of total particles.
Compare with the validation results in the following:
http://trappe.oit.umn.edu/
https://doi.org/10.1021/acs.jced.9b00756
https://doi.org/10.1002/aic.15816
"""

import argparse
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import physical_constants
from pyfeasst import macrostate_distribution

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--dictionary_input', type=json.loads, help='dictionary with the following parameters for each particle: temperature(Kelvin), expected vapor density(kg/m^3), expected liquid density(kg/m^3), total number of particles, and an optional xyz file name to initialize vapor and liquid',
            default='{"/feasst/particle/ethane.txt":{"temp":[236,256],"expect_vapor_dens":[20.19, 35.4],"expect_liquid_dens":[469.4, 434.42],"num_particles":[1000, 1000],"xyz_vapor":["",""],"xyz_liquid":["",""]},\
                    "/feasst/particle/n-butane.txt":{"temp":[354],"expect_vapor_dens":[32.35],"expect_liquid_dens":[496.9],"num_particles":[512],"xyz_vapor":["",""],"xyz_liquid":["",""]}}')
    parser.add_argument('--liquid_cutoff', type=float, default=14, help='real space cutoff distance in the liquid')
    parser.add_argument('--vapor_cutoff_frac_pbc', type=float, default=0.4, help='fraction of initial vapor pbc to set vapor cutoff')
    parser.add_argument('--tpc', type=int, default=int(1e5), help='trials per cycle')
    parser.add_argument('--equilibration_cycles', type=int, default=int(1e1),
                        help='number of cycles for equilibration')
    parser.add_argument('--production_cycles', type=int, default=int(1e2),
                        help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
    parser.add_argument('--run_type', '-r', type=int, default=0,
                        help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--seed', type=int, default=-1,
                        help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
    parser.add_argument('--scratch', type=str, default=None,
                        help='Optionally write scheduled job to scratch/logname/jobid.')
    parser.add_argument('--node', type=int, default=0, help='node ID')
    parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
    parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

    # Convert arguments into a parameter dictionary, and add argument-dependent parameters.
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    params = vars(args)
    params['script'] = __file__
    params['prefix'] = 'trappe'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['procs_per_sim'] = 1
    params['num_nodes'] = 1
    params['procs_per_node'] = 0
    params['particles'] = list()
    params['temperatures'] = list()
    params['vapor_pbcs'] = list()
    params['liquid_pbcs'] = list()
    params['num_vapors'] = list()
    params['num_liquids'] = list()
    params['xyz_vapors'] = list()
    params['xyz_liquids'] = list()
    params['molecular_weight'] = list()
    params['num_sitess'] = list()
    params['expect_vapor_dens'] = list()
    params['expect_liquid_dens'] = list()
    params['vapor_cutoffs'] = list()
    for _, part in enumerate(params['dictionary_input']):
        params['procs_per_node'] += len(params['dictionary_input'][part]['temp'])
        for index, temp in enumerate(params['dictionary_input'][part]['temp']):
            params['particles'].append(part)
            params['temperatures'].append(temp)
            params['num_sitess'].append(fstio.num_sites_in_fstprt(part, params['feasst_install']))
            if "ethane" in part:
                params['molecular_weight'].append(30.07)
            elif "n-butane" in part:
                params['molecular_weight'].append(58.12)
            else:
                print("input molecular weight for", part)
                assert False
            frac_vapor = 0.15
            num_part = params['dictionary_input'][part]['num_particles'][index]
            params['num_vapors'].append(int(num_part*frac_vapor))
            params['num_liquids'].append(num_part - params['num_vapors'][-1])
            dens_conv = density_convert(params['molecular_weight'][-1])
            params['expect_vapor_dens'].append(params['dictionary_input'][part]['expect_vapor_dens'][index])
            params['vapor_pbcs'].append((params['num_vapors'][-1]/params['expect_vapor_dens'][-1]*dens_conv)**(1./3.))
            params['expect_liquid_dens'].append(params['dictionary_input'][part]['expect_liquid_dens'][index])
            params['liquid_pbcs'].append((params['num_liquids'][-1]/params['expect_liquid_dens'][-1]*dens_conv)**(1./3.))
            params['vapor_cutoffs'].append(params['vapor_pbcs'][-1]*params['vapor_cutoff_frac_pbc'])
            params['xyz_vapors'].append(params['dictionary_input'][part]['xyz_vapor'][index])
            params['xyz_liquids'].append(params['dictionary_input'][part]['xyz_liquid'][index])
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    params['hours_terminate'] = 0.95*params['hours_terminate'] - 0.05 # terminate FEASST before SLURM
    params['hours_terminate'] *= params['procs_per_node'] # real time -> cpu time
    params['hours_checkpoint'] *= params['procs_per_node']
    #params['temperatures'] = np.linspace(params['temperature_lower'],params['temperature_upper'], num=params['num_sims']).tolist()
    params['mu_init']=10
    params['equil'] = params['equilibration_cycles']*params['tpc']
    params['double_equil'] = 2*params['equil']
    params['dccb_cut'] = 3.5 # cutoff for dual-cut cb in liquid
    params['num_dccb'] = 4   # number of cb steps per site in liquid
    return params, args

def density_convert(molecular_weight):
    """ multiply to convert N/V units of molecules/A^3 to kg/m^3 """
    na = physical_constants.AvogadroConstant().value()
    return 1./na*molecular_weight/1e3*1e30

def sim_node_dependent_params(params):
    """ Set parameters that depend upon the sim or node here. """
    sim = params['sim']
    params['beta'] = 1./(params['temperatures'][sim]*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
    params['vapor_pbc'] = params['vapor_pbcs'][sim]
    params['liquid_pbc'] = params['liquid_pbcs'][sim]
    params['num_vapor'] = params['num_vapors'][sim]
    params['num_liquid'] = params['num_liquids'][sim]
    params['xyz_vapor'] = params['xyz_vapors'][sim]
    params['xyz_liquid'] = params['xyz_liquids'][sim]
    params['fstprt'] = params['particles'][sim]
    params['num_sites'] = params['num_sitess'][sim]
    params['last_site'] = params['num_sites'] - 1
    params['vapor_cutoff'] = params['vapor_cutoffs'][sim]
    filepre = filename="{}{:03d}".format(params['prefix'], sim)
    fstio.write_linear_grow_file(filename=filepre+"_vapor_grow_canonical.txt", num_sites=params['num_sites'], gce=0, conf="vapor", ref="noixn", num_steps=1, particle_type='trappe')
    if params['xyz_vapor'] == '':
        fstio.write_linear_grow_file(filename=filepre+"_vapor_grow_add.txt", num_sites=params['num_sites'], gce=2, conf="vapor", ref="noixn", num_steps=1, particle_type='trappe')
    fstio.write_linear_grow_file(filename=filepre+"_liquid_grow_canonical.txt", num_sites=params['num_sites'], gce=0, conf="liquid", ref="dccb", num_steps=4, particle_type='trappe')
    if params['xyz_liquid'] == '':
        fstio.write_linear_grow_file(filename=filepre+"_liquid_grow_add.txt", num_sites=params['num_sites'], gce=2, conf="liquid", ref="dccb", num_steps=4, particle_type='trappe')
    fstio.write_linear_grow_file(filename=filepre+"_grow_gibbs.txt", num_sites=params['num_sites'], gce=3, conf="vapor", conf2="liquid", ref="dccb", num_steps=4, particle_type='trappe')

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed={seed}
For [config]:[cutoff]:[xyz]:[pbc]=vapor:{vapor_cutoff}:?{xyz_vapor}:?{vapor_pbc},liquid:{liquid_cutoff}:?{xyz_liquid}:?{liquid_pbc}
    Let [Config]=Configuration name=[config] particle_type=trappe:{fstprt} cutoff=[cutoff]
    If defined=?[xyz]
        [Config] xyz_file=[xyz]
    Else
        [Config] cubic_side_length=[pbc]
    Endif
    Potential Model=LennardJones config=[config]
    Potential Model=LennardJones VisitModel=VisitModelIntra intra_cut=3 config=[config]
    Potential VisitModel=LongRangeCorrections config=[config]
    RefPotential ref=noixn VisitModel=DontVisitModel config=[config]
EndFor
# Initialize dual-cut configuration bias reference potential in the liquid but not in the vapor
RefPotential ref=dccb config=vapor VisitModel=DontVisitModel
RefPotential ref=dccb config=liquid Model=LennardJones VisitModel=VisitModelCell cutoff={dccb_cut} min_length={dccb_cut}
ThermoParams beta={beta} chemical_potential={mu_init}
Metropolis
For [config]:[param]=vapor:30,liquid:1
    TrialTranslate weight=0.5 tunable_param=[param] config=[config]
    TrialGrowFile grow_file={prefix}{sim:03d}_[config]_grow_canonical.txt
EndFor
For [config]:[param]=vapor:180,liquid:0.4
    For [last_site]=0,{last_site}
        TrialParticlePivot weight=0.25 particle_type=trappe tunable_param=[param] pivot_site=[last_site] config=[config]
    EndFor
EndFor
CheckEnergy trials_per_update={tpc} decimal_places=4
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# gcmc initialization
Let [write]=trials_per_write={tpc} output_file={prefix}{sim:03d}
Log [write]_fill.csv
Tune
For [config]:[xyz]:[num]=vapor:?{xyz_vapor}:?{num_vapor},liquid:?{xyz_liquid}:?{num_liquid}
    Movie [write]_[config]_fill.xyz config=[config]
    If undefined=[xyz]
        TrialGrowFile grow_file={prefix}{sim:03d}_[config]_grow_add.txt
        Run until_num_particles=[num] config=[config]
        Remove name_contains=add
    EndIf
    Remove name=Movie
EndFor
Remove name=Tune,Log

# gibbs equilibration cycles: equilibrate, estimate density, adjust, repeat
# start a very long run GibbsInitialize completes once targets are reached
Metropolis trials_per_cycle=1e9 cycles_to_complete=1e9
GibbsInitialize updates_density_equil={equil} updates_per_adjust={double_equil}
TrialGrowFile grow_file={prefix}{sim:03d}_grow_gibbs.txt
TrialGibbsVolumeTransfer weight=0.006 tunable_param=3000 ref=noixn print_num_accepted=true configs=vapor,liquid
# a new tune is required when new Trials are introduced
# decrease trials per due to infrequency of volume transfer attempts
Tune trials_per_tune=20
Log [write]_eq.csv
For [config]=vapor,liquid
    Movie [write]_[config]_eq.xyz config=[config]
EndFor
ProfileCPU [write]_eq_profile.csv
# decrease trials per due to infrequency of volume transfer attempts
Run until=complete
Remove name=GibbsInitialize,Tune,Log,Movie,Movie,ProfileCPU

# gibbs ensemble production
Metropolis trials_per_cycle={tpc} cycles_to_complete={production_cycles}
Log [write].csv
For [analyze]:[file]=Density:_dens.csv,Movie:.xyz,Energy:_en.csv,Volume:_vol.csv,ProfileCPU:_profile.csv,CPUTime:_cpu.csv
    For [config]=vapor,liquid
        [analyze] [write]_[config][file] config=[config]
    EndFor
EndFor
GhostTrialVolume trials_per_update=1e3 [write]_pressure.csv
Run until=complete
""".format(**params))

def post_process(params):
    z_factor = 13
    for sim in range(params['num_sims']):
        part = params['particles'][sim]
        temp = params['temperatures'][sim]
        vapor_density = pd.read_csv("{}{:03d}_vapor_dens.csv".format(params['prefix'], sim))
        dens_conv = density_convert(params['molecular_weight'][sim])
        vapor_density['average'] *= dens_conv
        vapor_density['stdev'] *= dens_conv
        vapor_density['block_stdev'] *= dens_conv
        vapor_density['diff'] = np.abs(vapor_density['average']-params['expect_vapor_dens'][sim])
        vapor_density['tol'] = np.sqrt(vapor_density['block_stdev']**2+(2**2))
        print(part, temp, 'K vapor density', vapor_density)
        diverged = vapor_density[vapor_density['diff'] > z_factor*vapor_density['tol']]
        if len(diverged) > 0:
            print(diverged)
        assert len(diverged) == 0
        liquid_density = pd.read_csv("{}{:03d}_vapor_dens.csv".format(params['prefix'], sim))
        liquid_density['average'] *= dens_conv
        liquid_density['stdev'] *= dens_conv
        liquid_density['block_stdev'] *= dens_conv
        liquid_density['diff'] = np.abs(liquid_density['average']-params['expect_liquid_dens'][sim])
        liquid_density['tol'] = np.sqrt(liquid_density['block_stdev']**2+(2**2))
        print(part, temp, 'K liquid_density', liquid_density)
        diverged = liquid_density[liquid_density['diff'] > z_factor*liquid_density['tol']]
        if len(diverged) > 0:
            print(diverged)
        assert len(diverged) == 0
        if temp == 354 and "n-butane" in part:
            na = physical_constants.AvogadroConstant().value()
            pres_conv = 1e33/na # convert from kJ/mol/A^3 to Pa (J/m^3)
            pressure = pd.read_csv("{}{:03d}_pressure.csv".format(params['prefix'], sim))
            pressure['average'] *= pres_conv
            pressure['block_stdev'] *= pres_conv
            #pressure['diff'] = np.abs(pressure['average']-1.1976E+06)
            #pressure['tol'] = np.sqrt(pressure['block_stdev']**2+(3.6212E+02)**2)
            print(part, temp, 'K pressure', pressure)

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
