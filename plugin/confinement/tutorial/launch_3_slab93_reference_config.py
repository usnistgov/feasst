"""
Compute the energy for a single particle confined in a slab to test and debug the interaction model.
"""

import argparse
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/particle/lj.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--xy_side_length', type=float, default=9,
                    help='periodic boundary length in x and y (confined z)')
PARSER.add_argument('--z_side_length', type=float, default=6,
                    help='boundary length in z')
PARSER.add_argument('--z_particle_position', type=float, default=2.5,
                    help='boundary length in z')
PARSER.add_argument('--hours_terminate', type=float, default=0.2, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--seed', type=int, default=-1,
                    help='Random number generator seed. If -1, assign random seed to each sim.')
PARSER.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
PARSER.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
PARSER.add_argument('--scratch', type=str, default=None,
                    help='Optionally write scheduled job to scratch/logname/jobid.')
PARSER.add_argument('--node', type=int, default=0, help='node ID')
PARSER.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
PARSER.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

# Convert arguments into a parameter dictionary, and add argument-dependent parameters.
ARGS, UNKNOWN_ARGS = PARSER.parse_known_args()
assert len(UNKNOWN_ARGS) == 0, 'An unknown argument was included: '+str(UNKNOWN_ARGS)
PARAMS = vars(ARGS)
PARAMS['script'] = __file__
PARAMS['prefix'] = 'slab_ref_config'
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.95*PARAMS['hours_terminate'] - 0.05 # terminate FEASST before SLURM
PARAMS['hours_terminate'] *= PARAMS['procs_per_node'] # real time -> cpu time
PARAMS['num_sims'] = PARAMS['num_nodes']
PARAMS['procs_per_sim'] = PARAMS['procs_per_node']
PARAMS['half_z_side_length'] = PARAMS['z_side_length']/2.
with open(PARAMS['prefix']+'_shape_file.txt', 'w') as file1:
    file1.write("""Slab dimension 2 bound0 -{half_z_side_length} bound1 {half_z_side_length}""".format(**PARAMS))
with open(PARAMS['prefix']+'.xyz', 'w') as file1:
    file1.write("""1
-1 8 8 8 0 0 0
0 0 0 {z_particle_position}""".format(**PARAMS))

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
MonteCarlo
Configuration side_length0 {xy_side_length} side_length1 {xy_side_length} side_length2 {z_side_length} periodic2 false particle_type0 {fstprt} add_particles_of_type0 1 xyz_file {prefix}.xyz
Potential Model ModelLJShape shape_file {prefix}_shape_file.txt alpha 9 wall_epsilon 10 wall_sigma 2
Potential Model ModelLJShape shape_file {prefix}_shape_file.txt alpha 3 wall_epsilon -10 wall_sigma 2
Potential Model LennardJones
ThermoParams beta 1 chemical_potential 1
Metropolis
Log output_file {prefix}.csv max_precision true clear_file true
Run num_trials 1
""".format(**params))

def post_process(params):
    """ Compare the energy with the expected value """
    df = pd.read_csv(params['prefix']+'.csv')
    print(df)
    epsilon_mix = np.sqrt(10.*1.)
    sigma_mix = (2.+1.)/2.
    rcut = 3
    def en(r, epsilon, sigma):
        return epsilon*((sigma/r)**9 - (sigma/r)**3)
    shift = en(rcut, epsilon_mix, sigma_mix)
    print('shift', shift)
    assert np.abs(shift - -0.389108383966031) < 1e-15
    en_expected = en(0.5, epsilon_mix, sigma_mix) - shift
    print('en_expected', en_expected, 'diff', en_expected - df['energy'][0])
    assert np.abs(en_expected - df['energy'][0]) < 1e-15

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
