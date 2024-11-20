"""
Flat-histogram simulation of bulk CO2 to compare with the ZIF8 adsorption and obtain pressure.
Compare with https://doi.org/10.1021/jp400480q .
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import physical_constants

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fluid', type=str, default='/feasst/particle/co2.fstprt',
                    help='FEASST particle definition of a fluid particle / adsorbate.')
PARSER.add_argument('--cutoff', type=float, default=15, help='site-site cutoff distance in Angstroms')
PARSER.add_argument('--temperature', type=float, default=303, help='temperature in Kelvin')
PARSER.add_argument('--mu', type=float, default=-16, help='chemical potential')
PARSER.add_argument('--mu_init', type=float, default=-1, help='initial chemical potential')
PARSER.add_argument('--beta_init', type=float, default=0.1, help='initial beta to fill adsorbent')
PARSER.add_argument('--max_particles', type=int, default=395, help='maximum number of particles')
PARSER.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
PARSER.add_argument('--min_sweeps', type=int, default=1,
                    help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
PARSER.add_argument('--cubic_side_length', type=float, default=30,
                    help='cubic periodic boundary length')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e5),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=1e0,
                    help='number of iterations for equilibration')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.5, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=32, help='number of processors')
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
PARAMS['script'] = __file__
PARAMS['prefix'] = 'co2_'
PARAMS['alpha'] = 5.6/PARAMS['cubic_side_length']
PARAMS['beta'] = 1./(PARAMS['temperature']*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.95*PARAMS['hours_terminate'] - 0.05 # terminate FEASST before SLURM
PARAMS['hours_terminate'] *= PARAMS['procs_per_node'] # real time -> cpu time
PARAMS['hours_checkpoint'] *= PARAMS['procs_per_node']
PARAMS['num_sims'] = PARAMS['num_nodes']
PARAMS['procs_per_sim'] = PARAMS['procs_per_node']

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_checkpoint} ln_prob_file {prefix}{node}_lnpi.txt min_window_size -1
WindowExponential maximum {max_particles} min0 {min_particles} num {procs_per_node} overlap 0 alpha 2.2 min_size 2
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fluid} cutoff {cutoff}
Potential VisitModel Ewald alpha {alpha} kmax_squared 27
Potential Model ModelTwoBodyFactory model0 LennardJonesForceShift model1 ChargeScreened erfc_table_size 2e4 VisitModel VisitModelCutoffOuter
Potential Model ChargeScreenedIntra VisitModel VisitModelBond
Potential Model ChargeSelf
ThermoParams beta {beta_init} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 30 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 20 tunable_param 0.5 tunable_target_acceptance 0.25 particle_type 0
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-4

# gcmc initialization and nvt equilibration
TrialAdd weight 50.0 particle_type 0
Log trials_per_write {trials_per_iteration} output_file {prefix}{node}s[sim_index]_eq.txt
Tune
Run until_num_particles [soft_macro_min]
Remove name TrialAdd
ThermoParams beta {beta} chemical_potential {mu}
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Run until_criteria_complete true
Remove name0 Tune name1 Log

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
  Bias TransitionMatrix min_sweeps {min_sweeps}
TrialTransfer weight 50 particle_type 0
Log trials_per_write {trials_per_iteration} output_file {prefix}{node}s[sim_index].txt
Movie trials_per_write {trials_per_iteration} output_file {prefix}{node}s[sim_index]_eq.xyz stop_after_iteration 1
Movie trials_per_write {trials_per_iteration} output_file {prefix}{node}s[sim_index].xyz start_after_iteration 1
Tune trials_per_write {trials_per_iteration} output_file {prefix}{node}s[sim_index]_tune.txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per_iteration} output_file {prefix}{node}s[sim_index]_en.txt multistate true start_after_iteration 1
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {trials_per_iteration} output_file {prefix}{node}s[sim_index]_crit.txt
""".format(**params))

def post_process(params):
    assert True # placeholder

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
