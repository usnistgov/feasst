"""
This is a simple test of a henry coefficient calculation.
A particle is inserted in a box.
Half the box leads to overlap (e.g., a hard wall potential).
The <e^-(beta U)> term should be 1/2.
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
PARSER.add_argument('--beta', type=float, default=1., help='inverse temperature')
PARSER.add_argument('--mu', type=float, default=-1, help='chemical potential')
PARSER.add_argument('--cubic_side_length', type=float, default=9,
                    help='cubic periodic boundary length')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e5),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.02, help='hours per checkpoint')
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
PARAMS['prefix'] = 'henry'
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.95*PARAMS['hours_terminate'] - 0.05 # terminate FEASST before SLURM
PARAMS['hours_terminate'] *= PARAMS['procs_per_node'] # real time -> cpu time
PARAMS['hours_checkpoint'] *= PARAMS['procs_per_node']
PARAMS['num_sims'] = PARAMS['num_nodes']
PARAMS['procs_per_sim'] = PARAMS['procs_per_node']
with open(PARAMS['prefix']+'_shape_file.txt', 'w') as file1:
    file1.write("""HalfSpace dimension 0 intersection -0.5 direction 1""".format(**PARAMS))

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt}
Potential Model HardSphere
Potential Model ModelHardShape shape_file {prefix}_shape_file.txt
ThermoParams beta {beta} chemical_potential {mu}
AlwaysReject
TrialAdd particle_type 0 new_only true
HenryCoefficient trials_per_write {trials_per_iteration} file_name {prefix}.csv write_precision 12 num_beta_taylor 4
Run num_trials 1e6
""".format(**params))

def post_process(params):
    df = pd.read_csv(params['prefix'] + '.csv', comment='#')
    print(df)
    assert np.abs(df['average'][0] - 0.5) < 3*df['block_stdev'][0]
    with open(params['prefix']+'.csv') as f:
        firstline = f.readline().rstrip()
        henry=eval(firstline[1:])
        print(henry)
        assert np.abs(henry['beta_taylor'][0] - 0.5) < 3*df['block_stdev'][0]
        assert np.abs(henry['beta_taylor'][1]) < 1e-6

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
