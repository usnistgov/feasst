"""
Example Mayer-sampling (https://doi.org/10.1103/PhysRevLett.92.220601) simulation of a square well with a hard sphere reference.
Compare with Eq. 6 of https://doi.org/10.1063/1.1569473

Note that the post_process compares the simulated with the expected result using reduced which is the B2 divided by the B2 of the hard sphere of unit diameter.

An important aspect of computing second osmotic virial coefficients using Mayer-sampling is that you are calculating the ratio of the second virial of the full potential to the reference.
Thus, the reference potential needs to be something that you know already.
Typically, a hard sphere is a good choice.
In order to ensure that your reference potential is a single hard sphere, especially if your full potential has multiple sites, you will want to set model parameters in reference potential (sigma, cutoff) to zero.
Then you can set the sigma and cutoff of one particle site to the hard sphere size of interest.
Configuration sets parameters in a specific order, so that the setting of all parameters happens first, followed by the setting of particular parameters (see Configuration).

As an excercise, consider modifying this tutorial to compute the B2 of LennardJones with the following steps:
- Compare with a known result, such as temperature at which the B2 is zero: T_Boyle=3.417928023 from https://doi.org/10.1016/S0378-4371(00)00362-9
- Set beta = 1/T_Boyle
- Set cutoff to 1/2 the Domain side length (e.g., 5e9 for 1e10 box)
- Replace the line "Potential Model SquareWell" with "Potential Model LennardJones"
- Set b2reduced_analytical to zero.
- Compare with the lj tutorial.
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
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/atom.fstprt',
                        help='FEASST particle definition')
    parser.add_argument('--reference_sigma', type=float, default=1,
                        help='reference potential is a hard sphere unit diameter which is also the size of the inner hard sphere in the square well.')
    parser.add_argument('--cutoff', type=float, default=3,
                        help='the square well attractive interaction cutoff distance between centers')
    parser.add_argument('--beta', type=float, default=1./2., help='the inverse temperature')
    parser.add_argument('--tpc', type=int, default=int(1e4), help='trials per cycle')
    parser.add_argument('--equilibration_cycles', type=int, default=int(1e1),
                        help='number of cycles for equilibration')
    parser.add_argument('--production_cycles', type=int, default=int(1e2),
                        help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=0.1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=0.1, help='hours until termination')
    parser.add_argument('--procs_per_node', type=int, default=32, help='number of processors')
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
    params['prefix'] = 'sqw'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    return params, args

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length 1e10 periodic0 false periodic1 false periodic2 false \
  particle_type0 {fstprt} add_particles_of_type0 2 \
  group0 first first_particle_index 0 \
  cutoff {cutoff}
Potential Model SquareWell
RefPotential Model HardSphere sigma 0 sigma0 {reference_sigma} cutoff 0 cutoff0 {reference_sigma}
ThermoParams beta {beta}
MayerSampling trials_per_cycle {tpc} cycles_to_complete {equilibration_cycles}
TrialTranslate new_only true reference_index 0 tunable_param 1 group first
#TrialRotate new_only true reference_index 0 tunable_param 40
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# tune trial parameters
CriteriaWriter trials_per_write {tpc} output_file {prefix}{sim}_b2_eq.txt
Log trials_per_write {tpc} output_file {prefix}{sim}_eq.txt
Movie trials_per_write {tpc} output_file {prefix}{sim}_eq.xyz
Tune
Run until complete
Remove name0 Tune name1 CriteriaWriter name2 Log name3 Movie

# production
CriteriaWriter trials_per_write {tpc} output_file {prefix}{sim}_b2.txt
Log trials_per_write {tpc} output_file {prefix}{sim}.txt
Movie trials_per_write {tpc} output_file {prefix}{sim}.xyz
MayerSampling trials_per_cycle {tpc} cycles_to_complete {production_cycles}
Run until complete
""".format(**params))

def post_process(params):
    b2s = list()
    for sim in range(params['procs_per_node']):
        with open(params['prefix']+str(sim)+'_b2.txt') as f:
            firstline = f.readline().rstrip()
            b2=eval(firstline)
            #print(b2)
            b2s.append(b2['second_virial_ratio'])
    b2reduced_analytical = 1-(np.power(params['cutoff'], 3)-1)*(np.exp(params['beta'])-1)
    #b2hs = 2./3.*np.pi*params['reference_sigma']**3
    print('simulated', np.mean(b2s), 'std', np.std(b2s))
    print('expected', b2reduced_analytical)
    assert np.abs(np.mean(b2s) - b2reduced_analytical) < 5*np.std(b2s)

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
