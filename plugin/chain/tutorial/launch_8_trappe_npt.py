"""
Example NPT ensemble Monte Carlo simulation of TraPPE alkanes using FEASST.
Run multiple simulations, each on a single processor, using the argument "dictionary_input."
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
    parser.add_argument('--dictionary_input', type=json.loads, help='dictionary with the following parameters for each particle: temperature(Kelvin), initial density(kg/m^3), pressure(Pa), total number of particles, and an optional xyz file name to initialize the configuration',
#        default='{"/feasst/particle/ethane.txt":{"temp":[236],"initial_density":[20.19],"pressure":[1.1E+06],"num_particles":[512],"xyz":[""]}}')
#        default='{"/feasst/particle/n-butane.txt":{"temp":[354],"initial_density":[32.35],"pressure":[1.1976E+06],"num_particles":[512],"xyz":[""]}}')
        default='{"/feasst/particle/ethane.txt":{"temp":[236,256],"initial_density":[20.19, 35.4],"pressure":[1.1E+06,1.8E+06],"num_particles":[512, 512],"xyz":["",""]},\
                      "/feasst/particle/n-butane.txt":{"temp":[354],"initial_density":[32.35],"pressure":[1.1976E+06],"num_particles":[512],"xyz":[""]}}')
    parser.add_argument('--cutoff', type=float, default=14, help='real space cutoff distance')
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
    params['pressures'] = list()
    params['pbcs'] = list()
    params['nums'] = list()
    params['xyzs'] = list()
    params['molecular_weight'] = list()
    params['num_sitess'] = list()
    params['initial_density'] = list()
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
            num_part = params['dictionary_input'][part]['num_particles'][index]
            params['nums'].append(int(num_part))
            dens_conv = density_convert(params['molecular_weight'][-1])
            params['initial_density'].append(params['dictionary_input'][part]['initial_density'][index])
            params['pbcs'].append((params['nums'][-1]/params['initial_density'][-1]*dens_conv)**(1./3.))
            params['xyzs'].append(params['dictionary_input'][part]['xyz'][index])
            params['pressures'].append(params['dictionary_input'][part]['pressure'][index])
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    params['hours_terminate'] = 0.95*params['hours_terminate'] - 0.05 # terminate FEASST before SLURM
    params['hours_terminate'] *= params['procs_per_node'] # real time -> cpu time
    params['hours_checkpoint'] *= params['procs_per_node']
    #params['temperatures'] = np.linspace(params['temperature_lower'],params['temperature_upper'], num=params['num_sims']).tolist()
    params['mu_init']=10
    params['dccb_cut'] = 3.5 # cutoff for dual-cut cb
    params['num_dccb'] = 4   # number of cb steps per site
    return params, args

def density_convert(molecular_weight):
    """ multiply to convert N/V units of molecules/A^3 to kg/m^3 """
    na = physical_constants.AvogadroConstant().value()
    return 1./na*molecular_weight/1e3*1e30

def pressure_convert():
    """ multiply to convert from kJ/mol/A^3 to Pa (J/m^3) """
    na = physical_constants.AvogadroConstant().value()
    return 1e33/na

def sim_node_dependent_params(params):
    """ Set parameters that depend upon the sim or node here. """
    sim = params['sim']
    params['beta'] = 1./(params['temperatures'][sim]*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
    params['pbc'] = params['pbcs'][sim]
    params['num'] = params['nums'][sim]
    params['xyz'] = params['xyzs'][sim]
    params['pressure'] = params['pressures'][sim]/pressure_convert()
    params['fstprt'] = params['particles'][sim]
    params['num_sites'] = params['num_sitess'][sim]
    params['last_site'] = params['num_sites'] - 1
    if params['xyz'] == '':
        params['config'] = """cubic_side_length={pbc}""".format(**params)
        params['init'] = """TrialGrowFile grow_file={prefix}{sim:03d}_grow_add.txt
    Run until_num_particles={num}
    Remove name_contains=add""".format(**params)
    else:
        params['config'] = """xyz_file={xyz}""".format(**params)
        params['init'] = ''
    prepend="{}{:03d}".format(params['prefix'], sim)
    fstio.write_linear_grow_file(prepend+"_grow_canonical.txt", num_sites=params['num_sites'], gce=0, ref="noixn", num_steps=1, particle_type='trappe')
    if params['xyz'] == '':
        fstio.write_linear_grow_file(prepend+"_grow_add.txt", num_sites=params['num_sites'], gce=2, ref="noixn", num_steps=1, particle_type='trappe')

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed={seed}
Configuration {config} particle_type=trappe:{fstprt} cutoff={cutoff}
Potential Model=LennardJones VisitModel=VisitModelCell min_length=max_cutoff
Potential Model=LennardJones VisitModel=VisitModelIntra intra_cut=3
Potential VisitModel=LongRangeCorrections
RefPotential ref=noixn VisitModel=DontVisitModel
#RefPotential ref=dccb Model=LennardJones VisitModel=VisitModelCell cutoff={dccb_cut} min_length={dccb_cut}
ThermoParams beta={beta} chemical_potential={mu_init} pressure={pressure}
Metropolis
TrialTranslate weight=0.5 tunable_param=1
For [site]=0,{last_site}
    TrialParticlePivot weight=0.25 particle_type=trappe tunable_param=0.4 pivot_site=[site]
EndFor
TrialGrowFile grow_file={prefix}{sim:03d}_grow_canonical.txt
CheckEnergy trials_per_update={tpc} decimal_places=4
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# gcmc initialization
Let [write]=trials_per_write={tpc} output_file={prefix}{sim:03d}
Movie [write]_fill.xyz
Log [write]_fill.csv include_bonds=true
Tune
{init}
Remove name=Tune,Log,Movie

# nvt equilibration
Metropolis trials_per_cycle={tpc} cycles_to_complete={equilibration_cycles}
TrialVolume weight=0.005 tunable_param=0.2 tunable_target_acceptance=0.5
# decrease trials per due to infrequency of volume transfer attempts
Tune trials_per_tune=20
Log [write]_eq.csv
Movie [write]_eq.xyz
ProfileCPU [write]_eq_profile.csv
Run until=complete
Remove name=Tune,Log,Movie,ProfileCPU

# gibbs ensemble production
Metropolis trials_per_cycle={tpc} cycles_to_complete={production_cycles}
Log [write].csv
Movie [write].xyz
Energy [write]_en.csv
Density [write]_dens.csv
CPUTime [write]_cpu.csv
ProfileCPU [write]_profile.csv
Run until=complete
""".format(**params))

def post_process(params):
    z_factor = 13
    for sim in range(params['num_sims']):
        part = params['particles'][sim]
        temp = params['temperatures'][sim]
        density = pd.read_csv("{}{:03d}_dens.csv".format(params['prefix'], sim))
        dens_conv = density_convert(params['molecular_weight'][sim])
        density['average'] *= dens_conv
        density['stdev'] *= dens_conv
        density['block_stdev'] *= dens_conv
        density['diff'] = np.abs(density['average']-params['initial_density'][sim])
        density['tol'] = np.sqrt(density['block_stdev']**2+(2**2))
        print(part, temp, 'K density', density)
        diverged = density[density['diff'] > z_factor*density['tol']]
        if len(diverged) > 0:
            print(diverged)
        assert len(diverged) == 0

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
