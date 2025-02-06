"""
This is similar to the previous tutorial, except multiple MonteCarlo simulations
are run one after the other in the same text file.
"""

import argparse
import json
from pyfeasst import fstio

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/lj.fstprt',
                        help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=1./0.9, help='inverse temperature')
    parser.add_argument('--num_particles', type=int, default=500, help='number of particles')
    parser.add_argument('--pressures', type=json.loads, default='{"pressure":[8.9429E-04, 2.6485E-03, 4.3569E-03, 6.0193E-03, 7.6363E-03]}',
                        help='dictionary with a list of pressures to simulate')
    parser.add_argument('--initial_cubic_side_length', type=int, default=20, help='cubic periodic boundary length')
    parser.add_argument('--tpc', type=int, default=int(1e3), help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=int(1e1), help='number of cycles for equilibraiton')
    parser.add_argument('--production', type=int, default=int(1e1), help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
    parser.add_argument('--run_type', '-r', type=int, default=0,
                        help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--seed', type=int, default=-1,
                        help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
    parser.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
    parser.add_argument('--scratch', type=str, default=None,
                        help='Optionally write scheduled job to scratch/logname/jobid.')
    parser.add_argument('--queue_flags', type=str, default="", help='extra flags for queue (e.g., for slurm, "-p queue")')
    parser.add_argument('--node', type=int, default=0, help='node ID')
    parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
    parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

    # Convert arguments into a parameter dictionary, and add argument-dependent parameters.
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    params = vars(args)
    params['script'] = __file__
    params['prefix'] = 'mlj'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_node'] = 1
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    return params, args

def write_feasst_script(params, script_file):
    script = ""
    for sim, pressure in enumerate(params['pressures']['pressure']):
        print('pressure', pressure)
        params['pressure'] = pressure
        params['sim'] = sim
        script += one_sim(params) + '\n'
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write(script)

def one_sim(params):
    """ Write fst script for a single simulations """
    return """MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {initial_cubic_side_length} particle_type0 {fstprt}
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential -1
Metropolis
TrialTranslate tunable_param 2 tunable_target_acceptance 0.2
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}
CheckEnergy trials_per_update {tpc} decimal_places 4

# gcmc initialization
TrialAdd particle_type 0
Log trials_per_write {tpc} output_file {prefix}{sim}_init.csv
Tune
Run until_num_particles {num_particles}
Remove name0 TrialAdd name1 Log name2 Tune

# npt equilibration
ThermoParams beta {beta} pressure {pressure}
Metropolis trials_per_cycle {tpc} cycles_to_complete {equilibration}
TrialVolume weight 0.005 tunable_param 0.2 tunable_target_acceptance 0.5
Tune    trials_per_tune 20
Log     trials_per_write {tpc} output_file {prefix}{sim}_eq.csv
Movie   trials_per_write {tpc} output_file {prefix}{sim}_eq.xyz
Density trials_per_write {tpc} output_file {prefix}{sim}_density_eq.csv
Run until complete
Remove name0 Tune name1 Log name2 Movie name3 Density

# npt production
Metropolis trials_per_cycle {tpc} cycles_to_complete {production}
Log        trials_per_write {tpc} output_file {prefix}{sim}.csv
Movie      trials_per_write {tpc} output_file {prefix}{sim}.xyz
Energy     trials_per_write {tpc} output_file {prefix}{sim}_en.csv
Density    trials_per_write {tpc} output_file {prefix}{sim}_density.csv
Volume     trials_per_write {tpc} output_file {prefix}{sim}_volume.csv
ProfileCPU trials_per_write {tpc} output_file {prefix}{sim}_cpu.csv
Run until complete
""".format(**params)

def post_process(params):
    """ Plot energy and compare with https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    ens = np.zeros(shape=(params['num_sims'], 2))
    rhos = np.zeros(shape=(params['num_sims'], 2))
    for sim in range(params['num_sims']):
        log = pd.read_csv(params['prefix']+str(sim)+'.csv')
        assert int(log['num_particles_of_type0'][0]) == params['num_particles']
        energy = pd.read_csv(params['prefix']+str(sim)+'_en.csv')
        ens[sim] = np.array([energy['average'][0],
                             energy['block_stdev'][0]])/params['num_particles']
        density = pd.read_csv(params['prefix']+str(sim)+'_density.csv')
        print('density', density)
        rhos[sim] = np.array([density['average'][0],
                              density['block_stdev'][0]])
        print('rhos[sim]', rhos[sim])
    # data from https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
    #rhos_srsw = [0.001, 0.003, 0.005, 0.007, 0.009]
    print('rhos', rhos)
    #ens_srsw = [-9.9165E-03, -2.9787E-02]
    ens_srsw = [-9.9165E-03, -2.9787E-02, -4.9771E-02, -6.9805E-02, -8.9936E-02]
    en_stds_srsw = [1.89E-05, 3.21E-05]

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
