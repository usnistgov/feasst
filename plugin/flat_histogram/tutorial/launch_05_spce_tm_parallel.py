"""
Flat-histogram simulation of SPC/E water in the grand canonical ensemble.
The default temperature of 525 K is above the critical point.
Dual-cut configurational bias is used for insertions and deletions with a reference potential of
hard spheres on the oxygen sigma.
The results are are compared with the SRSW https://doi.org/10.18434/T4M88Q
(which use CODATA2010 physical constants).
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import feasstio
from pyfeasst import physical_constants

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default=os.path.expanduser('~')+'/feasst/build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/forcefield/spce.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--temperature', type=float, default=525, help='temperature in Kelvin')
PARSER.add_argument('--beta_mu', type=float, default=-8.14, help='beta times chemical potential')
PARSER.add_argument('--mu_init', type=float, default=10, help='initial chemical potential')
PARSER.add_argument('--max_particles', type=int, default=265, help='maximum number of particles')
PARSER.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
PARSER.add_argument('--min_sweeps', type=int, default=5,
                    help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
PARSER.add_argument('--cubic_box_length', type=float, default=20,
                    help='cubic periodic boundary length')
PARSER.add_argument('--dccb_cut', type=float, default=0.9*3.165,
                    help='dual-cut configurational bias cutoff')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e6),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=0,
                    help='number of iterations for equilibraiton')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.02, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=0.2, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=32, help='number of processors')
PARSER.add_argument('--prefix', type=str, default='spce', help='prefix for all output file names')
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
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.95*PARAMS['hours_terminate'] - 0.05 # terminate FEASST before SLURM
PARAMS['hours_terminate'] *= PARAMS['procs_per_node'] # real time -> cpu time
PARAMS['hours_checkpoint'] *= PARAMS['procs_per_node']
PARAMS['num_sims'] = PARAMS['num_nodes']
PARAMS['procs_per_sim'] = PARAMS['procs_per_node']
PARAMS['alpha'] = 5.6/PARAMS['cubic_box_length']
PARAMS['beta'] = 1./(PARAMS['temperature']*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
PARAMS['mu'] = PARAMS['beta_mu']/PARAMS['beta']
PARAMS['dccb_cut'] = PARAMS['cubic_box_length']/int(PARAMS['cubic_box_length']/PARAMS['dccb_cut']) # maximize inside box
PARAMS['place_h'] = 'angle true mobile_site 2 anchor_site 0 anchor_site2 1 reference_index 0'
if 'tip4p' in PARAMS['fstprt']:
    PARAMS['place_h'] = 'branch true mobile_site 2 mobile_site2 3 anchor_site 0 anchor_site2 1 reference_index 0'

# write TrialGrowFile for SPCE
with open('spce_grow.txt', 'w') as f:
    f.write("""TrialGrowFile

particle_type 0 weight 2 transfer true site 0 num_steps 10 reference_index 0
bond true mobile_site 1 anchor_site 0 reference_index 0
{place_h}
""".format(**PARAMS))

def write_feasst_script(params, file_name):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(file_name, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_checkpoint} ln_prob_file {prefix}n{node}_lnpi.txt min_window_size -1
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 1.5 min_size 5
Checkpoint file_name {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt} group0 oxygen oxygen_site_type 0
Potential VisitModel Ewald alpha {alpha} kmax_squared 38
Potential Model ModelTwoBodyFactory model0 LennardJones model1 ChargeScreened erfc_table_size 2e4 VisitModel VisitModelCutoffOuter
RefPotential Model HardSphere group oxygen cutoff {dccb_cut} VisitModel VisitModelCell min_length {dccb_cut} cell_group oxygen
Potential Model ChargeScreenedIntra VisitModel VisitModelBond
Potential Model ChargeSelf
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 0.5 particle_type 0 tunable_param 0.5 tunable_target_acceptance 0.25
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-4

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Log trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_eq.txt
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
Log trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index].txt
Movie trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_eq.xyz stop_after_iteration 1
Movie trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index].xyz start_after_iteration 1
Tune trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_tune.txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_en.txt multistate true start_after_iteration 1
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_crit.txt
""".format(**params))

def post_process(params):
    """ Skip the following checks if temperature is not 525 K """
    if np.abs(params['temperature'] - 525) > 1e-5:
        return
    lnpi = pd.read_csv(params['prefix']+'n0_lnpi.txt')
    lnpi = pd.concat([lnpi, pd.read_csv(params['feasst_install']+'../plugin/flat_histogram/test/data/stat_spce_525.csv')], axis=1)
    lnpi['deltalnPI'] = lnpi.lnPI - lnpi.lnPI.shift(1)
    diverged = lnpi[lnpi.deltalnPI-lnpi.delta_ln_prob > 6*lnpi.delta_ln_prob_stdev]
    print(diverged)
    assert len(diverged) == 0
    plt.plot(lnpi['state'], lnpi['ln_prob'], label='FEASST')
    plt.plot(lnpi['N'], lnpi['lnPI'], linestyle='dashed', label='SRSW')
    plt.xlabel('number of particles', fontsize=16)
    plt.ylabel('ln probability', fontsize=16)
    plt.legend(fontsize=16)
    #plt.savefig(params['prefix']+'_lnpi.png', bbox_inches='tight', transparent='True')

if __name__ == '__main__':
    feasstio.run_simulations(params=PARAMS,
                             sim_node_dependent_params=None,
                             write_feasst_script=write_feasst_script,
                             post_process=post_process,
                             queue_function=feasstio.slurm_single_node,
                             args=ARGS)
