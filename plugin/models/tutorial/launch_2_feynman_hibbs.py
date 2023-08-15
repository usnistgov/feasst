"""
This tutorial shows how to build your own LJ-like potentials, specifically, Feynman-Hibbs.
Note that the table potential may be faster because it avoids the sqrt operation.
"""

import os
import subprocess
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import feasstio

# Parse arguments from command line or change their default values.
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--feasst_install', type=str, default=os.path.expanduser('~')+'/feasst/build/',
    help='FEASST install directory (e.g., the path to build)')
parser.add_argument('--fstprt', type=str, default='/feasst/forcefield/lj.fstprt',
    help='FEASST particle definition')
parser.add_argument('--beta', type=float, default=1./0.9, help='inverse temperature')
parser.add_argument('--num_particles', type=int, default=50, help='number of particles')
parser.add_argument('--density', type=float, default=0.001, help='density')
parser.add_argument('--trials_per_iteration', type=int, default=int(1e5),
    help='like cycles, but not necessary num_particles')
parser.add_argument('--equilibration_iterations', type=int, default=int(1e1),
    help='number of iterations for equilibraiton')
parser.add_argument('--production_iterations', type=int, default=int(1e3),
    help='number of iterations for production')
parser.add_argument('--hours_checkpoint', type=float, default=0.1, help='hours per checkpoint')
parser.add_argument('--hours_terminate', type=float, default=0.1, help='hours until termination')
parser.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
parser.add_argument('--prefix', type=str, default='fh', help='prefix for all output file names')
parser.add_argument('--run_type', '-r', type=int, default=0,
    help='0: run, 1: submit to queue, 2: post-process')
parser.add_argument('--seed', type=int, default=-1,
    help='Random number generator seed. If -1, assign random seed to each sim.')
parser.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
parser.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
parser.add_argument('--scratch', type=str, default=None,
    help='Optionally write scheduled job to scratch/logname/jobid.')
PARSER.add_argument('--queue_flags', type=str, default="", help='extra flags for queue (e.g., for slurm, "-p queue")')
parser.add_argument('--node', type=int, default=0, help='node ID')
parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

# Convert arguments into a parameter dictionary, and add argument-dependent parameters.
args, unknown_args = parser.parse_known_args()
assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
params = vars(args)
epsilon = 1. #4*eps*h^2/24mkbT/sigma^2
params['s_14'] = 132*epsilon
params['s_8'] = -30*epsilon
params['cubic_box_length'] = (params['num_particles']/params['density'])**(1./3.)
params['script'] = __file__
params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
params['hours_terminate'] = 0.95*params['hours_terminate'] - 0.05 # terminate before queue
params['procs_per_sim'] = 1
params['num_sims'] = params['num_nodes']*params['procs_per_node']

def write_feasst_script(params, file_name):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(file_name, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt}
Potential Model ModelTwoBodyFactory \
  model0 LennardJones \
  model1 TwoBodyAlpha alpha0 14 s0 {s_14} alpha1 8 s1 {s_8} \
  VisitModel VisitModelCell min_length 3
ThermoParams beta {beta} chemical_potential -1
Metropolis
TrialTranslate tunable_param 2 tunable_target_acceptance 0.2
Checkpoint file_name {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# grand canonical ensemble initalization
TrialAdd particle_type 0
Run until_num_particles {num_particles}
RemoveTrial name TrialAdd

# canonical ensemble equilibration
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Tune
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-8
Log trials_per_write {trials_per_iteration} file_name {prefix}{sim}_eq.txt
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log

# canonical ensemble production
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {production_iterations}
Log trials_per_write {trials_per_iteration} file_name {prefix}{sim}.txt
Movie trials_per_write {trials_per_iteration} file_name {prefix}{sim}.xyz
Energy trials_per_write {trials_per_iteration} file_name {prefix}{sim}_en.txt
CPUTime trials_per_write {trials_per_iteration} file_name {prefix}{sim}_cpu.txt
Run until_criteria_complete true
""".format(**params))

def post_process(params):
    # test not implemented
    assert True

if __name__ == '__main__':
    feasstio.run_simulations(params=params,
                             sim_node_dependent_params=None,
                             write_feasst_script=write_feasst_script,
                             post_process=post_process,
                             queue_function=feasstio.slurm_single_node,
                             args=args)
