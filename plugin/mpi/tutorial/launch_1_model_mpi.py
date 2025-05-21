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

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/lj.txt',
                        help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=1/0.9, help='1 / kB / T')
    parser.add_argument('--density', type=float, default=1e-3, help='number density')
    parser.add_argument('--num_particles', type=int, default=256, help='number of particles')
    parser.add_argument('--tpc', type=int, default=int(1e5), help='trials per cycle')
    parser.add_argument('--equilibration_cycles', type=int, default=int(1e1),
                        help='number of cycles for equilibraiton')
    parser.add_argument('--production_cycles', type=int, default=int(1e1),
                        help='number of cycles for production')
    parser.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
    parser.add_argument('--seed', type=int, default=-1,
                        help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--run_type', '-r', type=int, default=0,
                        help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--hours_terminate', type=float, default=0.1, help='hours until termination')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
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
    params['prefix'] = 'lj'
    params['script'] = __file__
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_sim'] = 1
    params['num_sims'] = params['procs_per_node']
    params['cubic_side_length'] = np.power(params['num_particles']/params['density'], 1./3.)
    return params, args

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt}
Potential Model ModelMPI
#Potential Model ModelMPI VisitModel VisitModelCell
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
    parameters, arguments = parse()
    if args.run_type == 1:
        client(params)
    elif args.run_type == 0:
        params['sim'] = 0
        write_feasst_script(params, params['prefix']+"0_run.txt")
        cmd="""mpirun -np 1 {feasst_install}bin/fst < {prefix}0_run.txt : -np 1 python {script} -r 1""".format(**params)
        print('cmd:', cmd)
        #subprocess.call(cmd, shell=True, executable='/bin/bash')
        from subprocess import Popen
        process = Popen(cmd, shell=True)
        post_process(params)
    elif args.run_type == 2:
        post_process(params)
#    fstio.run_simulations(params=params,
#                          queue_function=fstio.slurm_single_node,
#                          args=args,
#                          write_feasst_script=write_feasst_script,
#                          client=client,
#                          post_process=post_process)
