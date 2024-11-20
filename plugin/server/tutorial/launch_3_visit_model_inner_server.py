"""
Define site-site energy from a python client (ModelServer) for generic rigid bodies in 3D.
"""

import subprocess
import multiprocessing
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import argparse
import random
import socket
from pyfeasst import fstio

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/plugin/aniso/particle/aniso_tabular.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--beta', type=float, default=1/0.9, help='1 / kB / T')
PARSER.add_argument('--port', type=int, default=50005, help='server client interface port')
PARSER.add_argument('--density', type=float, default=1e-3, help='number density')
PARSER.add_argument('--num_particles', type=int, default=256, help='number of particles')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e5),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=int(1e1),
                    help='number of iterations for equilibraiton')
PARSER.add_argument('--production_iterations', type=int, default=int(1e1),
                    help='number of iterations for production')
PARSER.add_argument('--buffer_size', type=int, default=1000, help='server client interface port')
PARSER.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
PARSER.add_argument('--seed', type=int, default=-1,
                    help='Random number generator seed. If -1, assign random seed to each sim.')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--hours_terminate', type=float, default=0.1, help='hours until termination')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.1, help='hours per checkpoint')
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
PARAMS['prefix'] = 'aniso'
PARAMS['script'] = __file__
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.99*PARAMS['hours_terminate'] - 0.0333 # terminate before queue
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['procs_per_node']
PARAMS['cubic_side_length'] = np.power(PARAMS['num_particles']/PARAMS['density'], 1./3.)

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt}
Potential Model TwoBodyTable VisitModelInner VisitModelInnerServer port {port} server_site0 0
ThermoParams beta {beta} chemical_potential -1
Metropolis
TrialTranslate tunable_param 2 tunable_target_acceptance 0.2
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}
TrialAdd particle_type 0
Run until_num_particles {num_particles}
Remove name TrialAdd
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Tune
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-8
Log trials_per_write {trials_per_iteration} output_file {prefix}{sim}_eq.csv
Run until_criteria_complete true
Remove name0 Tune name1 Log
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {production_iterations}
Log trials_per_write {trials_per_iteration} output_file {prefix}{sim}.csv
Movie trials_per_write {trials_per_iteration} output_file {prefix}{sim}.xyz
Energy trials_per_write {trials_per_iteration} output_file {prefix}{sim}_en.csv
CPUTime trials_per_write {trials_per_iteration} output_file {prefix}{sim}_cpu.txt
Run until_criteria_complete true
""".format(**params))

def en_lj(r2):
    return 4*(1./r2**6 -1./r2**3)

def client(params):
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.connect(("localhost", params['port']+params['sim']))
    terminate = False
    while not terminate:
        message = sock.recv(1000)
        if len(message) == 0:
            terminate = True
        else:
            decoded_string = message.decode("utf-8")
            vals = decoded_string.split(',')
            r2 = float(vals[0])
            s1 = float(vals[1])
            s2 = float(vals[2])
            e1 = float(vals[3])
            e2 = float(vals[4])
            e3 = float(vals[5])
            type1 = int(vals[6])
            type2 = int(vals[7])
            #print(r2, s1, s2, e1, e2, e3, type1, type2)
            sock.send(str(en_lj(r2)).encode('utf-8'))
    sock.close()

def post_process(params):
    assert True # place holder

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS,
                          write_feasst_script=write_feasst_script,
                          client=client,
                          post_process=post_process)
