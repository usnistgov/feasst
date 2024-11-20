"""
Conduct ghost insertions during a two-macrostate flat-histogram simulation to compare with flat histogram transition probabilities.
"""

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
PARSER.add_argument('--beta', type=float, default=1./1.5, help='inverse temperature')
PARSER.add_argument('--mu', type=float, default=-2.352321, help='chemical potential')
PARSER.add_argument('--mu_init', type=float, default=10, help='initial chemical potential')
PARSER.add_argument('--max_particles', type=int, default=21, help='maximum number of particles')
PARSER.add_argument('--min_particles', type=int, default=20, help='minimum number of particles')
PARSER.add_argument('--min_sweeps', type=int, default=1e3,
                    help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
PARSER.add_argument('--cubic_side_length', type=float, default=8,
                    help='cubic periodic boundary length')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e6),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=0,
                    help='number of iterations for equilibration')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.02, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=0.15, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--seed', type=int, default=-1,
                    help='Random number generator seed. If -1, assign random seed to each sim.')
PARSER.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
PARSER.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
PARSER.add_argument('--scratch', type=str, default=None,
                    help='Optionally write scheduled job to scratch/logname/jobid.')
PARSER.add_argument('--node', type=int, default=0, help='node ID')
PARSER.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
PARSER.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

# Convert arguments into a parameter dictionary, and add argument-dependent parameters.
ARGS, UNKNOWN_ARGS = PARSER.parse_known_args()
assert len(UNKNOWN_ARGS) == 0, 'An unknown argument was included: '+str(UNKNOWN_ARGS)
PARAMS = vars(ARGS)
PARAMS['prefix'] = 'ghost'
PARAMS['script'] = __file__
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.95*PARAMS['hours_terminate'] - 0.05 # terminate FEASST before SLURM
PARAMS['hours_terminate'] *= PARAMS['procs_per_node'] # real time -> cpu time
PARAMS['hours_checkpoint'] *= PARAMS['procs_per_node']
PARAMS['num_sims'] = PARAMS['num_nodes']
PARAMS['procs_per_sim'] = PARAMS['procs_per_node']

# write TrialGrowFile
with open(PARAMS['prefix']+'_ghost_insertion.txt', 'w') as f:
    f.write("""TrialGrowFile

particle_type 0 add true site 0 num_steps 10
""")

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_checkpoint} ln_prob_file {prefix}n{node}_lnpi.txt min_window_size -1
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 2.25 min_size 2
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt}
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
CheckEnergy trials_per_update {trials_per_iteration} decimal_places 4

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Log trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_eq.txt
Tune
Run until_num_particles [soft_macro_min]
Remove name TrialAdd
ThermoParams beta {beta} chemical_potential {mu}
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Run until_criteria_complete true
Remove name0 Tune name1 Log

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} min_flatness 25 collect_flatness 20 min_collect_sweeps 1 widom true
TrialTransfer weight 2 particle_type 0
Log trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index].txt
#To print trajectories for each macrostate in separate files, add the following arguments to the "Movie" lines below: multistate true multistate_aggregate false
Movie trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_eq.xyz stop_after_iteration 1
Movie trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index].xyz start_after_iteration 1
Tune trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_tune.txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_en.txt multistate true start_after_iteration 1
GhostTrialGrow trials_per_update 10 trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_ghost.txt grow_file {prefix}_ghost_insertion.txt multistate true
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_crit.txt
""".format(**params))

def post_process(params):
    """ An analysis is not yet implemented """
    assert True

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
