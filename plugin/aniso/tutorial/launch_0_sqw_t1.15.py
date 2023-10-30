"""
This tutorial is similar to flat histogram tutorial 4, but for a square well
with anisotropic "texture" from random noise.
This is more of a mechanical test of the anisotropic potentials reproducing
an isotropic potential:
https://www.nist.gov/mml/csd/chemical-informatics-group/sat-tmmc-liquid-vapor-coexistence-properties-square-well-fluid
"""

import argparse
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import macrostate_distribution

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/plugin/aniso/particle/aniso_tabular.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--cutoff', type=float, default=1.5, help='square well cutoff')
PARSER.add_argument('--dccb_cut', type=float, default=1, help='DCCB cutoff')
PARSER.add_argument('--num_orientations_per_pi', type=int, default=3,
                    help='number of orientations per 90 degrees in each table angle')
PARSER.add_argument('--num_z', type=int, default=2, help='number of distances in table')
PARSER.add_argument('--table_file', type=str, default='dat.txt', help='table file name')
PARSER.add_argument('--beta', type=float, default=1./1.15, help='inverse temperature')
PARSER.add_argument('--mu', type=float, default=-3, help='chemical potential')
PARSER.add_argument('--mu_init', type=float, default=10, help='initial chemical potential')
PARSER.add_argument('--max_particles', type=int, default=475, help='maximum number of particles')
PARSER.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
PARSER.add_argument('--min_sweeps', type=int, default=20,
                    help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
PARSER.add_argument('--cubic_side_length', type=float, default=9,
                    help='cubic periodic boundary length')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e6),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=0,
                    help='number of iterations for equilibration')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.02, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=0.2, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=32, help='number of processors')
PARSER.add_argument('--prefix', type=str, default='sqw', help='prefix for all output file names')
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
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.95*PARAMS['hours_terminate'] - 0.05 # terminate FEASST before SLURM
PARAMS['hours_terminate'] *= PARAMS['procs_per_node'] # real time -> cpu time
PARAMS['hours_checkpoint'] *= PARAMS['procs_per_node']
PARAMS['num_sims'] = PARAMS['num_nodes']
PARAMS['procs_per_sim'] = PARAMS['procs_per_node']

# generate the tabular potential
def generate_table(num_orientations_per_pi, num_z, table_file, gamma=1., delta=0.5):
    dt = np.pi/num_orientations_per_pi
    assert num_z > 1
    dz = 1./(num_z - 1)
    num_orientations = 0
    num_elements = 0
    fl = open(table_file, 'w')
    fl.write('site_types 1 0\n' +
             'num_orientations_per_pi ' + str(num_orientations_per_pi) + '\n' +
             'gamma ' + str(gamma) + '\n' +
             'delta ' + str(delta) + '\n' +
             'num_z ' + str(num_z) + '\n' +
             'smoothing_distance -1\n')
    #for s1 in np.arange(0, 2*np.pi + dt/2, dt): #theta
    for s1 in np.arange(0, np.pi + dt/2, dt): # if i==j, avoid x < 0 for i/j swap symmetry
        for s2 in np.arange(0, np.pi + dt/2, dt): #phi
            for e1 in np.arange(-np.pi, np.pi + dt/2, dt):
                for e2 in np.arange(0, np.pi + dt/2, dt):
                    for e3 in np.arange(-np.pi, np.pi + dt/2, dt):
                        num_orientations += 1
                        # hard code square well for testing
                        rh = 1.
                        fl.write(str(rh) + " ")
                        for z in np.arange(0, 1. + dz/2, dz):
                            num_elements += 1
                            # hard code square well for testing
                            energy = -1. + 0.0001*random.uniform(-1., 1.)
                            fl.write(str(energy) + " ")
                        fl.write("\n")
    fl.close()
    print('num_orientations', num_orientations)
    print('num_elements', num_elements)

generate_table(num_orientations_per_pi=ARGS.num_orientations_per_pi,
               num_z=ARGS.num_z,
               table_file=ARGS.table_file)

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_checkpoint} ln_prob_file {prefix}n{node}_lnpi.txt min_window_size -1
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 2.25 min_size 5
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt} cutoff {cutoff}
Potential Model TwoBodyTable VisitModelInner VisitModelInnerTable table_file dat.txt
RefPotential Model HardSphere cutoff {dccb_cut} VisitModel VisitModelCell min_length {dccb_cut}
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
TrialRotate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-4

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Log trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_eq.txt
Tune
Run until_num_particles [soft_macro_min]
RemoveTrial name TrialAdd
ThermoParams beta {beta} chemical_potential {mu}
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} min_flatness 25 collect_flatness 20 min_collect_sweeps 1
TrialTransfer weight 2 particle_type 0
Log trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index].txt
Movie trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_eq.xyz stop_after_iteration 1
Movie trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index].xyz start_after_iteration 1
Tune trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_tune.txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_en.txt multistate true start_after_iteration 1
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_crit.txt
""".format(**params))

def post_process(params):
    lnpi = macrostate_distribution.MacrostateDistribution(file_name=params['prefix']+'n0_lnpi.txt')
    rw = lnpi.equilibrium()
    self.assertAlmostEqual(params['beta']*params['mu'] + rw, -3.194, delta=1e-2)
    #lnpi.plot(show=True)
    vap, liq = lnpi.split()
    assert np.abs(vap.average_macrostate()/params['cubic_side_length']**3 - 9.723E-02) < 1e-3
    assert np.abs(liq.average_macrostate()/params['cubic_side_length']**3 - 5.384E-01) < 1e-3

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
