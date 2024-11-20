"""
Tabulate a potential with multiple site types.
At the default temperature, the trimer model, also demonstrated in the flat_histogram tutorial 8,
forms micelles which have difficulty sampling in this canonical ensemble simulation,
especially when there are no cluster moves (see https://doi.org/10.1063/1.4918557).
Regardless, the formation of micelle-like structures demonstrates the tabular potential may be
appropriately accounting for the different site types.
Note that if more than one particle is present, the site types indices add to the ones of the
previous type (e.g., 2 particles with two site types each would have site types 0 and 1 in
particle type 0 and site types 2 and 3 in particle type 1.
"""

import subprocess
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/particle/trimer.fstprt',
    help='FEASST particle definition')
PARSER.add_argument('--beta', type=float, default=1./0.2, help='inverse temperature')
PARSER.add_argument('--num_particles', type=int, default=500, help='number of particles')
PARSER.add_argument('--rho_lower', type=float, default=9e-2, help='lowest number density')
PARSER.add_argument('--rho_upper', type=float, default=9e-2, help='highest number density')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e5),
    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=int(1e1),
    help='number of iterations for equilibration')
PARSER.add_argument('--production_iterations', type=int, default=int(1e2),
    help='number of iterations for production')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.1, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=0.1, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
PARSER.add_argument('--table_file', type=str, default='trimer_table.txt', help='table file name')
PARSER.add_argument('--num_z', type=int, default=int(1e3), help='number of table elements')
PARSER.add_argument('--inner', type=float, default=0.75, help='As described in TablePotential')
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
PARAMS['prefix'] = 'trimer'
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.99*PARAMS['hours_terminate'] - 0.0333 # terminate before queue
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']
PARAMS['rhos'] = np.linspace(PARAMS['rho_lower'], PARAMS['rho_upper'], num=PARAMS['num_sims'])
PARAMS['cubic_side_lengths'] = np.power(PARAMS['num_particles']/PARAMS['rhos'], 1./3.).tolist()
PARAMS['rhos'] = PARAMS['rhos'].tolist()
PARAMS['gamma'] = -2 # as described in TablePotential
PARAMS['cutoff'] = 3
PARAMS['rwca'] = 2**(1/6.)

def sim_node_dependent_params(params):
    params['cubic_side_length'] = params['cubic_side_lengths'][params['sim']]

def user_potential(distance, site1, site2):
    en_lj = 4*(distance**-12 - distance**-6)
    if site1 == 0 and site2 == 0:
        return en_lj
    else:
        return en_lj + 1

def generate_table(PARAMS):
    assert PARAMS['num_z'] > 1
    dz = 1./(PARAMS['num_z'] - 1)
    with open(PARAMS['table_file'], 'w') as file1:
        file1.write("""site_types 2 0 1""")
        for site1 in [0, 1]:
            for site2 in [0, 1]:
                if site1 <= site2:
                    if site1 == 0 and site2 == 0:
                        cutoff = PARAMS['cutoff']
                    else:
                        cutoff = PARAMS['rwca']
                    rcg = cutoff**PARAMS['gamma']
                    rhg = PARAMS['inner']**PARAMS['gamma']
                    file1.write("""\ninner {inner}\nnum_z {num_z}\n""".format(**PARAMS))
                    #file1.write("""site_types 1 0\ngamma {gamma}\ninner {inner}\nnum_z {num_z}\n""".format(**PARAMS))
                    for z in np.arange(0, 1 + dz/2, dz):
                        if z == 0:
                            distance = PARAMS['inner']
                        else:
                            distance = (z*(rcg - rhg) + rhg)**(1./PARAMS['gamma'])
                        en = user_potential(distance, site1, site2)
                        #print('distance', distance, 'en', en)
                        file1.write(str(en) + " ")
generate_table(PARAMS)

# write fst script to run a simulation
def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt} cutoff {cutoff} cutoff0_1 {rwca} cutoff1_1 {rwca}
Potential Model TablePotential table_file {table_file}
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential -1
Metropolis
TrialTranslate tunable_param 2 tunable_target_acceptance 0.2
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# grand canonical ensemble initalization
TrialAdd particle_type 0
Run until_num_particles {num_particles}
Remove name TrialAdd

# canonical ensemble equilibration
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Tune
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-8
Log trials_per_write {trials_per_iteration} output_file {prefix}{sim}_eq.txt
Run until_criteria_complete true
Remove name0 Tune name1 Log

# canonical ensemble production
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {production_iterations}
Log trials_per_write {trials_per_iteration} output_file {prefix}{sim}.txt
Movie trials_per_write {trials_per_iteration} output_file {prefix}{sim}.xyz
Energy trials_per_write {trials_per_iteration} output_file {prefix}{sim}_en.txt
CPUTime trials_per_write {trials_per_iteration} output_file {prefix}{sim}_cpu.txt
Run until_criteria_complete true
""".format(**params))

def post_process(params):
    """ placeholder """
    assert True

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
