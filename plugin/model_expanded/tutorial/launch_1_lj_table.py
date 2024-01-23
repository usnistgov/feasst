"""
This tutorial is essentially a copy of /path/to/feasst/tutorial/launch.py,
except that the LJ potential is tabulated with different sigma values,
and an expanded ensemble is considered in sigma.

Normalize distances by the cubic side length, V=1^3.
"""

import subprocess
import argparse
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
PARSER.add_argument('--beta', type=float, default=1./0.9, help='inverse temperature')
PARSER.add_argument('--num_particles', type=int, default=500, help='number of particles')
PARSER.add_argument('--rho', type=float, default=1e-3, help='lowest number density')
PARSER.add_argument('--eps_lower', type=float, default=0.5, help='lowest number density')
PARSER.add_argument('--eps_upper', type=float, default=1, help='highest number density')
PARSER.add_argument('--eps_num', type=int, default=6, help='number of rhoes from lower to upper')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e5),
    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=int(1e0),
    help='number of iterations for equilibration')
PARSER.add_argument('--min_flatness', type=int, default=25, help='number of WL flatness')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.1, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=0.1, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
PARSER.add_argument('--num_z', type=int, default=int(1e3), help='number of table elements')
PARSER.add_argument('--inner', type=float, default=0.75, help='As described in TablePotential')
PARSER.add_argument('--cutoff', type=float, default=3, help='potential cutoff distance')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--seed', type=int, default=-1,
    help='Random number generator seed. If -1, assign random seed to each sim.')
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
PARAMS['prefix'] = 'lj'
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.99*PARAMS['hours_terminate'] - 0.0333 # terminate before queue
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']
PARAMS['epsilons'] = np.linspace(PARAMS['eps_lower'], PARAMS['eps_upper'], num=PARAMS['eps_num']).tolist()
PARAMS['cubic_side_length'] = np.power(PARAMS['num_particles']/PARAMS['rho'], 1./3.)
PARAMS['max'] = len(PARAMS['epsilons']) - 1
PARAMS['gamma'] = -2 # as described in TablePotential

def user_potential(distance, epsilon, sigma=1):
    return 4*epsilon*((sigma/distance)**12 - (sigma/distance)**6)

def generate_table(params):
    assert params['num_z'] > 1
    dz = 1./(params['num_z'] - 1)
    rhg = params['inner']**params['gamma']
    rcg = params['cutoff']**params['gamma']
    with open(params['table_file'], 'w') as file1:
        file1.write("""site_types 1 0\ninner {inner}\nnum_z {num_z}\n""".format(**params))
        #file1.write("""site_types 1 0\ngamma {gamma}\ninner {inner}\nnum_z {num_z}\n""".format(**params))
        for z in np.arange(0, 1 + dz/2, dz):
            if z == 0:
                distance = params['inner']
            else:
                distance = (z*(rcg - rhg) + rhg)**(1./params['gamma'])
            en = user_potential(distance=distance, epsilon=params['epsilon'])
            #print('distance', distance, 'en', en)
            file1.write(str(en) + " ")

with open(PARAMS['prefix']+'_models.txt', 'w') as file1:
    file1.write('ModelTwoBodyFactory\n\n')
    for index, epsilon in enumerate(PARAMS['epsilons']):
        PARAMS['epsilon'] = epsilon
        PARAMS['table_file'] = """{prefix}e{epsilon}.txt""".format(**PARAMS)
        generate_table(PARAMS)
        file1.write("""TablePotential table_file {table_file}\n""".format(**PARAMS))

# write fst script to run a simulation
def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt}
Potential Model ModelExpanded model_file {prefix}_models.txt VisitModel VisitModelCell min_length max_cutoff
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential -1
Metropolis
TrialTranslate tunable_param 2 tunable_target_acceptance 0.2
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# grand canonical ensemble initalization
TrialAdd particle_type 0
Run until_num_particles {num_particles}
RemoveTrial name TrialAdd

# canonical ensemble equilibration
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Tune
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-8
Log trials_per_write {trials_per_iteration} output_file {prefix}{sim}_eq.txt
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log

# canonical ensemble production
FlatHistogram Macrostate MacrostateModel width 1 max {max} min 0 Bias WangLandau min_flatness {min_flatness}
TrialModel
CriteriaUpdater trials_per_update {trials_per_iteration}
CriteriaWriter trials_per_write {trials_per_iteration} output_file {prefix}{sim}_crit.txt
Log trials_per_write {trials_per_iteration} output_file {prefix}{sim}.txt
Movie trials_per_write {trials_per_iteration} output_file {prefix}{sim}.xyz
Energy trials_per_write {trials_per_iteration} output_file {prefix}{sim}_en.txt multistate true
CPUTime trials_per_write {trials_per_iteration} output_file {prefix}{sim}_cpu.txt
Run until_criteria_complete true
""".format(**params))

def post_process(params):
    """ Plot energy and compare with https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm """
    ens = pd.read_csv(params['prefix']+'0_en.txt')
    print(ens)
    print(ens['average']/params['num_particles'])
    for macro in range(len(ens)):
        # data from https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
        en_srsw = -9.9165E-03*params['epsilons'][macro]
        en_std_srsw = 1.89E-05*params['epsilons'][macro]
        print('macro', macro, 'srsw', en_srsw, en_std_srsw)
        diff = ens['average'][macro]/params['num_particles'] - en_srsw
        print('diff', diff)
        tol = 10*np.sqrt((ens['block_stdev'][macro]/params['num_particles'])**2 + en_std_srsw**2)
        print('tol', tol)
        assert np.abs(diff) < tol

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
