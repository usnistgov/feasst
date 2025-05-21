"""
Example Mayer-sampling (https://doi.org/10.1103/PhysRevLett.92.220601)
simulation of a Lennard-Jones particle while extrapolating in inverse
temperature, as described in https://doi.org/10.1063/1.5016165.
The reference potential is a hard sphere.
This reproduces Figure 2 of https://doi.org/10.1063/1.5016165.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import fstplot

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/atom.txt',
                        help='FEASST particle definition')
    parser.add_argument('--reference_sigma', type=float, default=1,
                        help='reference potential is a hard sphere unit diameter which is also the size of the inner hard sphere in the square well.')
    parser.add_argument('--beta', type=float, default=1., help='the inverse temperature')
    parser.add_argument('--num_beta_taylor', type=int, default=10, help='number of Taylor series coefficients')
    parser.add_argument('--show_plot', type=int, default=0, help='If != 0, show plot')
    parser.add_argument('--tpc', type=int, default=int(1e4), help='trials per cycle')
    parser.add_argument('--equilibration_cycles', type=int, default=int(1e1),
                        help='number of cycles for equilibration')
    parser.add_argument('--production_cycles', type=int, default=int(1e2),
                        help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
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
    params['prefix'] = 'lj'
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
  cutoff 5e9
Potential Model LennardJones
RefPotential Model HardSphere sigma 0 sigma0 {reference_sigma} cutoff 0 cutoff0 {reference_sigma}
ThermoParams beta {beta}
MayerSampling num_beta_taylor {num_beta_taylor} trials_per_cycle {tpc} cycles_to_complete {equilibration_cycles}
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
MayerSampling num_beta_taylor {num_beta_taylor} trials_per_cycle {tpc} cycles_to_complete {production_cycles}
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
            if sim == 0:
                df = pd.DataFrame(data={sim: b2['beta_taylor']})
            else:
                df[sim] = b2['beta_taylor']
    print('b2', np.mean(b2s), '+/-', np.std(b2s)/np.sqrt(len(b2s)))
    b2hs = 2./3.*np.pi*params['reference_sigma']**3
    coeffs = b2hs*df.mean(axis=1)
    print('Taylor series coefficients:', coeffs)
    coeffs_expected = [-5.322731267117734, -9.281005175555862, -2.810611771043566, -0.649558510117600]
    for index, co_exp in enumerate(coeffs_expected):
        if np.abs(coeffs[index] - co_exp) > 0.2:
            print(coeffs[index], co_exp)
            assert False
    if params['show_plot']:
        from scipy.interpolate import pade
        xrng = np.arange(0.01, 3.5, 0.1)
        beta0 = 1.
        orders = range(1, params['num_beta_taylor']+1)
        colors = fstplot.val2map(orders)
        for order, co in enumerate(coeffs):
            if order > 0:
                p_pade, q_pade = pade(list(coeffs[:(order+1)]), 1)
                deta = xrng - beta0
                plt.plot(xrng, p_pade(deta)/q_pade(deta), color=colors.to_rgba(order), marker='o')
        plt.xlim([0, 3.0])
        plt.ylim([-50, 10])
        fstplot.display(orders, label='order')
        plt.show()

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
