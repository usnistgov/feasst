"""
Example of running feasst as a server and interacting with a python client.
In this example, the simulation is terminated by checking elapsed time on the Python client.
"""

import subprocess
import multiprocessing
import time
import argparse
import random
import socket
from pyfeasst import fstio

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--beta', type=float, default=1/1.5, help='1 / kB / T')
    parser.add_argument('--port', type=int, default=54321, help='server client interface port')
    parser.add_argument('--buffer_size', type=int, default=1000, help='server client interface port')
    parser.add_argument('--procs_per_node', type=int, default=2, help='number of processors')
    parser.add_argument('--seed', type=int, default=-1,
                        help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--run_type', '-r', type=int, default=0,
                        help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--hours_terminate', type=float, default=0.04, help='hours until termination')
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
    params['prefix'] = 'listen'
    params['script'] = __file__
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_sim'] = 1
    params['num_sims'] = params['procs_per_node']
    return params, args

def client(params):
    start_time = time.time()
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.connect(("localhost", params['port']+params['sim']))
    feasst_commands = """MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length 8 particle_type0 /feasst/particle/lj.txt
Potential Model LennardJones
ThermoParams beta {beta} chemical_potential0 1
Metropolis
TrialTranslate
TrialAdd particle_type 0
Run until_num_particles 20
Remove name TrialAdd
Log trials_per_write 2e6 output_file {prefix}{sim}.csv
Run num_trials 2e6""".format(**params)
    for line in feasst_commands.split('\n'):
        sock.send(bytes(line, 'utf-8'))
        message = sock.recv(params['buffer_size'])
    print('elapsed', time.time() - start_time)
    print('terminate', params['hours_terminate']*60*60)
    while time.time() - start_time < params['hours_terminate']*60*60:
        sock.send(bytes('Run num_trials 1e5', 'utf-8'))
        message = sock.recv(params['buffer_size'])
        print('elapsed', params['sim'], time.time() - start_time)
    sock.close()

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          queue_function=fstio.slurm_single_node,
                          args=arguments,
                          client=client)
