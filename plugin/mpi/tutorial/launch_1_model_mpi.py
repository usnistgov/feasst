"""
Define site-site energy from a python client (ModelMPI).
"""

import subprocess
import multiprocessing
from mpi4py import MPI
import struct
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
PARSER.add_argument('--fstprt', type=str, default='/feasst/particle/lj.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--beta', type=float, default=1/0.9, help='1 / kB / T')
PARSER.add_argument('--density', type=float, default=1e-3, help='number density')
PARSER.add_argument('--num_particles', type=int, default=256, help='number of particles')
PARSER.add_argument('--tpc', type=int, default=int(1e5), help='trials per cycle')
PARSER.add_argument('--equilibration_cycles', type=int, default=int(1e1),
                    help='number of cycles for equilibraiton')
PARSER.add_argument('--production_cycles', type=int, default=int(1e1),
                    help='number of cycles for production')
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
PARAMS['prefix'] = 'lj'
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
Potential Model ModelMPI
#Potential Model ModelMPI VisitModel VisitModelCell min_length max_cutoff
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential -1
Metropolis
TrialTranslate tunable_param 2 tunable_target_acceptance 0.2
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}
TrialAdd particle_type 0
Run until_num_particles {num_particles}
Remove name TrialAdd
Metropolis trials_per_cycle {tpc} cycles_to_complete {equilibration_cycles}
Tune
CheckEnergy trials_per_update {tpc} tolerance 1e-8
Log trials_per_write {tpc} output_file {prefix}{sim}_eq.csv
Run until complete
Remove name0 Tune name1 Log
Metropolis trials_per_cycle {tpc} cycles_to_complete {production_cycles}
Log trials_per_write {tpc} output_file {prefix}{sim}.csv
Movie trials_per_write {tpc} output_file {prefix}{sim}.xyz
Energy trials_per_write {tpc} output_file {prefix}{sim}_en.csv
CPUTime trials_per_write {tpc} output_file {prefix}{sim}_cpu.txt
Run until complete
""".format(**params))

def en_lj(r2):
    return 4*(1./r2**6 -1./r2**3)

def client(params):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    print('rank', rank)
    r2_buffer = bytearray(b" " * 8)
    type1_buffer = bytearray(b" " * 4)
    type2_buffer = bytearray(b" " * 4)
    terminate = False
    while not terminate:
        comm.Recv([r2_buffer, MPI.FLOAT], source=0)
        r2 = struct.unpack('d', r2_buffer)[0]
#        print('r2', r2)
        if r2 > 0:
            comm.Recv([type1_buffer, MPI.INT], source=0)
            comm.Recv([type2_buffer, MPI.INT], source=0)
            #print('type1', struct.unpack('i', type1_buffer)[0])
            #print('type2', struct.unpack('i', type2_buffer)[0])
            en = en_lj(r2)
            comm.Send([bytearray(struct.pack('d', en)), MPI.FLOAT], dest=0)
        else:
            terminate = True

def post_process(params):
    """ Plot energy and compare with https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm """
    ens = np.zeros(shape=(params['num_sims'], 2))
    for sim in range(params['num_sims']):
        log = pd.read_csv(params['prefix']+str(sim)+'.csv')
        assert int(log['num_particles_of_type0'][0]) == params['num_particles']
        energy = pd.read_csv(params['prefix']+str(sim)+'_en.csv')
        ens[sim] = np.array([energy['average'][0],
                             energy['block_stdev'][0]])/params['num_particles']
    rhos_srsw = [0.001]
    ens_srsw = [-9.9165E-03]
    en_stds_srsw = [1.89E-05]
    print('ens', ens)
    plt.errorbar(rhos_srsw, ens_srsw, en_stds_srsw, fmt='+', label='SRSW')
    plt.errorbar(params['density'], ens[:, 0], ens[:, 1], fmt='x', label='FEASST')
    plt.xlabel(r'$\rho$', fontsize=16)
    plt.ylabel(r'$U/N$', fontsize=16)
    plt.legend(fontsize=16)
    #plt.savefig(params['prefix']+'_energy.png', bbox_inches='tight', transparent='True')
    if len(rhos_srsw) == params['num_sims']: # compare with srsw exactly
        for sim in range(params['num_sims']):
            diff = ens[sim][0] - ens_srsw[sim]
            assert np.abs(diff) < 10*np.sqrt(ens[sim][1]**2 + en_stds_srsw[sim]**2)

if __name__ == '__main__':
    if ARGS.run_type == 1:
        client(PARAMS)
    elif ARGS.run_type == 0:
        PARAMS['sim'] = 0
        write_feasst_script(PARAMS, PARAMS['prefix']+"0_run.txt")
        cmd="""mpirun -np 1 {feasst_install}bin/fst < {prefix}0_run.txt : -np 1 python {script} -r 1""".format(**PARAMS)
        print('cmd:', cmd)
        #subprocess.call(cmd, shell=True, executable='/bin/bash')
        from subprocess import Popen
        process = Popen(cmd, shell=True)
        post_process(PARAMS)
    elif ARGS.run_type == 2:
        post_process(PARAMS)
#    fstio.run_simulations(params=PARAMS,
#                          queue_function=fstio.slurm_single_node,
#                          args=ARGS,
#                          write_feasst_script=write_feasst_script,
#                          client=client,
#                          post_process=post_process)
