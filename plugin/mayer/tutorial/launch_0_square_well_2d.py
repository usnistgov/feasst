"""
Example Mayer-sampling (https://doi.org/10.1103/PhysRevLett.92.220601) simulation of a 2D square well with a hard circle reference.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/atom2d.txt',help='FEASST particle definition')
    parser.add_argument('--reference_sigma', type=float, default=1, help='reference potential is a hard sphere unit diameter which is also the size of the inner hard sphere in the square well.')
    parser.add_argument('--cutoff', type=float, default=1.5, help='the square well attractive interaction cutoff distance between centers')
    parser.add_argument('--beta', type=float, default=1./2., help='the inverse temperature')
    parser.add_argument('--tpc', type=int, default=int(1e4), help='trials per cycle')
    parser.add_argument('--equilibration_cycles', type=int, default=int(1e1),
                        help='number of cycles for equilibration')
    parser.add_argument('--production_cycles', type=int, default=int(1e2),
                        help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
    parser.add_argument('--num_jobs', type=int, default=3, help='Number of jobs in queue')
    parser.add_argument('--procs_per_job', type=int, default=1, help='number of processors')
    parser.add_argument('--run_type', '-r', type=int, default=0,
                        help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--seed', type=int, default=-1,
                        help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
    parser.add_argument('--scratch', type=str, default=None,
                        help='Optionally write scheduled job to scratch/logname/jobid.')
    parser.add_argument('--queue_flags', type=str, default="", help='extra flags for queue (e.g., for slurm, "-p queue")')
    parser.add_argument('--job', type=int, default=0, help='job ID')
    parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
    parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

    # Convert arguments into a parameter dictionary, and add argument-dependent parameters.
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    params = vars(args)
    params['script'] = __file__
    params['prefix'] = 'sqw2d'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_jobs']*params['procs_per_job']
    return params, args

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed={seed}
Configuration side_length=1e10,1e10 periodic=false,false \
    particle_type=fluid:{fstprt} add_num_fluid_particles=2 \
    group=first first_particle_index=0 \
    cutoff={cutoff}
Potential Model=SquareWell
RefPotential ref=hs Model=HardSphere sigma=0 sigmaA={reference_sigma} cutoff=0 cutoffA={reference_sigma}
ThermoParams beta={beta}
MayerSampling trials_per_cycle={tpc} cycles_to_complete={equilibration_cycles}
TrialTranslate new_only=true ref=hs tunable_param=1 group=first
#TrialRotate new_only=true ref=hs tunable_param=40
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# tune trial parameters
Let [write]=trials_per_write={tpc} output_file={prefix}{sim:03d}
Log [write]_eq.csv
Movie [write]_eq.xyz
CriteriaWriter [write]_b2_eq.csv
Tune
Run until=complete
Remove name=Tune,CriteriaWriter,Log,Movie

# production
Log [write].csv
Movie [write].xyz
CriteriaWriter [write]_b2.csv
MayerSampling trials_per_cycle={tpc} cycles_to_complete={production_cycles}
Run until=complete
""".format(**params))

def post_process(params):
    b2s = list()
    for sim in range(params['num_sims']):
        with open("{}{:03d}_b2.csv".format(params['prefix'], sim)) as f:
            firstline = f.readline().rstrip()
            b2=eval(firstline)
            #print(b2)
            b2s.append(b2['second_virial_ratio'])
    b2reduced_analytical = 1-(np.power(params['cutoff'], 2)-1)*(np.exp(params['beta'])-1)
    #b2hs = 2./3.*np.pi*params['reference_sigma']**3
    print('simulated', np.mean(b2s), 'std', np.std(b2s))
    print('expected', b2reduced_analytical)
    assert np.abs(np.mean(b2s) - b2reduced_analytical) < 5*np.std(b2s)

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_job_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_job,
                          args=arguments)
