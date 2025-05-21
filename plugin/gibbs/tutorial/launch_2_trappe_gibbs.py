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
    if params['xyz_vapor'] == '':
        params['vapor_config'] = """cubic_side_length {vapor_pbc}""".format(**params)
        params['init_vapor'] = """TrialGrowFile grow_file {prefix}{sim}_c0_grow_add.txt
    Run until_num_particles {num_vapor} configuration_index 0
    Remove name_contains add""".format(**params)
    else:
        params['vapor_config'] = """xyz_file {xyz_vapor}""".format(**params)
        params['init_vapor'] = ''
    if params['xyz_liquid'] == '':
        params['liquid_config'] = """cubic_side_length {liquid_pbc}""".format(**params)
        params['init_liquid'] = """TrialGrowFile grow_file {prefix}{sim}_c1_grow_add.txt
    Run until_num_particles {num_liquid} configuration_index 1
    Remove name_contains add""".format(**params)
    else:
        params['liquid_config'] = """xyz_file {xyz_liquid}""".format(**params)
        params['init_liquid'] = ''
    fstio.write_linear_grow_file(filename=params['prefix']+str(sim)+"_c0_grow_canonical.txt", num_sites=params['num_sites'], gce=0, conf=0, reference_index=0, num_steps=1)
    if params['xyz_vapor'] == '':
        fstio.write_linear_grow_file(filename=params['prefix']+str(sim)+"_c0_grow_add.txt", num_sites=params['num_sites'], gce=2, conf=0, reference_index=0, num_steps=1)
    fstio.write_linear_grow_file(filename=params['prefix']+str(sim)+"_c1_grow_canonical.txt", num_sites=params['num_sites'], gce=0, conf=1, reference_index=1, num_steps=4)
    if params['xyz_liquid'] == '':
        fstio.write_linear_grow_file(filename=params['prefix']+str(sim)+"_c1_grow_add.txt", num_sites=params['num_sites'], gce=2, conf=1, reference_index=1, num_steps=4)
    fstio.write_linear_grow_file(filename=params['prefix']+str(sim)+"_grow_gibbs.txt", num_sites=params['num_sites'], gce=3, conf=0, conf2=1, reference_index=1, num_steps=4)

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration {vapor_config} particle_type0 {fstprt} cutoff {vapor_cutoff}
Configuration {liquid_config} particle_type0 {fstprt} cutoff {liquid_cutoff}
CopyFollowingLines for_num_configurations 2
    Potential Model LennardJones
    Potential Model LennardJones VisitModel VisitModelIntra intra_cut 3
    Potential VisitModel LongRangeCorrections
    RefPotential reference_index 0 VisitModel DontVisitModel
EndCopy
# Initialize dual-cut configuration bias reference potential in the liquid but not in the vapor
RefPotential reference_index 1 configuration_index 0 VisitModel DontVisitModel
RefPotential reference_index 1 configuration_index 1 Model LennardJones VisitModel VisitModelCell cutoff {dccb_cut} min_length {dccb_cut}
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
CopyNextLine replace0 configuration_index with0 0 replace1 tunable_param with1 30
TrialTranslate weight 0.5 tunable_param 1 configuration_index 1
CopyNextLine replace0 configuration_index with0 0 replace1 tunable_param with1 180
TrialParticlePivot weight 0.25 particle_type 0 tunable_param 0.4 pivot_site 0 configuration_index 1
CopyNextLine replace0 configuration_index with0 0 replace1 tunable_param with1 180
TrialParticlePivot weight 0.25 particle_type 0 tunable_param 0.4 pivot_site {last_site} configuration_index 1
TrialGrowFile grow_file {prefix}{sim}_c0_grow_canonical.txt
TrialGrowFile grow_file {prefix}{sim}_c1_grow_canonical.txt
CheckEnergy trials_per_update {tpc} decimal_places 4
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# gcmc initialization and nvt equilibration
CopyNextLine replace0 configuration_index with0 0 replace1 output_file with1 {prefix}{sim}_c0_fill.xyz
Movie trials_per_write {tpc} output_file {prefix}{sim}_c1_fill.xyz configuration_index 1
Log trials_per_write {tpc} output_file {prefix}{sim}_fill.csv include_bonds true
Tune

# fill the first box
{init_vapor}

# fill the second box
{init_liquid}

Remove name0 Tune name1 Log name2 Movie name3 Movie

# gibbs equilibration cycles: equilibrate, estimate density, adjust, repeat
# start a very long run GibbsInitialize completes once targets are reached
Metropolis trials_per_cycle 1e9 cycles_to_complete 1e9
GibbsInitialize updates_density_equil {equil} updates_per_adjust {double_equil}
TrialGrowFile grow_file {prefix}{sim}_grow_gibbs.txt
TrialGibbsVolumeTransfer weight 0.006 tunable_param 3000 reference_index 0 print_num_accepted true
# a new tune is required when new Trials are introduced
# decrease trials per due to infrequency of volume transfer attempts
Tune trials_per_tune 20
CheckEnergy trials_per_update {tpc} decimal_places 8
Log trials_per_write {tpc} output_file {prefix}{sim}_eq.csv
CopyNextLine replace0 configuration_index with0 0 replace1 output_file with1 {prefix}{sim}_c0_eq.xyz
Movie trials_per_write {tpc} output_file {prefix}{sim}_c1_eq.xyz configuration_index 1
ProfileCPU trials_per_write {tpc} output_file {prefix}{sim}_eq_profile.csv
# decrease trials per due to infrequency of volume transfer attempts
Run until complete
Remove name0 GibbsInitialize name1 Tune name2 Log name3 Movie name4 Movie name5 ProfileCPU

# gibbs ensemble production
Metropolis trials_per_cycle {tpc} cycles_to_complete {production_cycles}
Log trials_per_write {tpc} output_file {prefix}{sim}.csv
CopyFollowingLines for_num_configurations 2 replace_with_index [config]
    Density trials_per_write {tpc} output_file {prefix}{sim}_c[config]_dens.csv
    Movie   trials_per_write {tpc} output_file {prefix}{sim}_c[config].xyz
    Energy  trials_per_write {tpc} output_file {prefix}{sim}_c[config]_en.csv
EndCopy
GhostTrialVolume trials_per_update 1e3 trials_per_write {tpc} output_file {prefix}{sim}_pressure.csv
ProfileCPU trials_per_write {tpc} output_file {prefix}{sim}_profile.csv
CPUTime    trials_per_write {tpc} output_file {prefix}{sim}_cpu.csv
Run until complete
""".format(**params))

def post_process(params):
    z_factor = 13
    for sim in range(params['num_sims']):
        part = params['particles'][sim]
        temp = params['temperatures'][sim]
        vapor_density = pd.read_csv(params['prefix']+str(sim)+"_c0_dens.csv")
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
        liquid_density = pd.read_csv(params['prefix']+str(sim)+"_c1_dens.csv")
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
            pressure = pd.read_csv(params['prefix']+str(sim)+"_pressure.csv")
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
