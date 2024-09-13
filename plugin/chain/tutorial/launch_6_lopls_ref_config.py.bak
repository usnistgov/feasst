"""
L-OPLS n-hexane reference configuration
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/plugin/chain/particle/n-hexane_lopls.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--cubic_side_length', type=float, default=203.4876,
                    help='cubic periodic boundary conditions')
PARSER.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--seed', type=int, default=-1,
                    help='Random number generator seed. If -1, assign random seed to each sim.')
PARSER.add_argument('--hours_terminate', type=float, default=1., help='hours until termination')
PARSER.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
PARSER.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
PARSER.add_argument('--scratch', type=str, default=None,
                    help='Optionally write scheduled job to scratch/logname/jobid.')
PARSER.add_argument('--queue_flags', type=str, default="", help='extra flags for queue (e.g., for slurm, "-p queue")')
PARSER.add_argument('--node', type=int, default=0, help='node ID')
PARSER.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
PARSER.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

# Convert arguments into a parameter dictionary, and add argument-dependent parameters.
ARGS, UNKNOWN_ARGS = PARSER.parse_known_args()
assert len(UNKNOWN_ARGS) == 0, 'An unknown argument was included: '+str(UNKNOWN_ARGS)
PARAMS = vars(ARGS)
PARAMS['script'] = __file__
PARAMS['prefix'] = 'lopls'
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.95*PARAMS['hours_terminate'] - 0.05 # terminate before queue
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = 1
PARAMS['alpha'] = 5.6/PARAMS['cubic_side_length']
assert False # this tutorial is incomplete
print('*****\nWARNING: ChargeScreenedIntra may not have the right VisitModel\n*****')

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt} add_particles_of_type0 1
Potential Model LennardJones
Potential Model LennardJones VisitModel VisitModelIntraMap exclude_bonds true exclude_angles true dihedral_weight 0.5
#Potential VisitModel Ewald alpha {alpha} kmax_squared 38
#Potential Model ModelTwoBodyFactory model0 LennardJones model1 ChargeScreened erfc_table_size 2e4
#Potential Model ChargeScreenedIntra VisitModel VisitModelBond
#Potential Model ChargeSelf
#Potential VisitModel LongRangeCorrections
ThermoParams beta 1 chemical_potential -1
Metropolis
Log output_file {prefix}{sim}.csv include_bonds true max_precision true clear_file true
Run num_trials 1
""".format(**params))

def post_process(params):
    df = pd.read_csv(params['prefix']+'0.csv')
    print(df)

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
