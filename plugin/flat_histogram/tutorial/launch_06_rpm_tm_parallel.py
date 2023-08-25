"""
Flat-histogram simulation of monovalent RPM in the grand canonical ensemble.
Dual-cut configurational bias is used for insertions and deletions.
The results are are compared with https://doi.org/10.1063/1.5123683
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import feasstio
from pyfeasst import physical_constants

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--plus', type=str, default='/feasst/plugin/charge/particle/rpm_plus.fstprt',
                    help='FEASST particle definition of the positive charge RPM')
PARSER.add_argument('--minus', type=str, default='/feasst/plugin/charge/particle/rpm_minus.fstprt',
                    help='FEASST particle definition of the negative charge RPM')
PARSER.add_argument('--beta', type=float, default=1./0.047899460618081, help='inverse temperature')
PARSER.add_argument('--beta_mu', type=float, default=-13.94, help='beta times chemical potential')
PARSER.add_argument('--max_particles', type=int, default=10, help='maximum number of particles')
PARSER.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
PARSER.add_argument('--min_sweeps', type=int, default=1e2,
                    help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
PARSER.add_argument('--cubic_side_length', type=float, default=12,
                    help='cubic periodic boundary length')
PARSER.add_argument('--dccb_cut', type=float, default=2**(1./6.),
                    help='dual-cut configurational bias cutoff')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e6),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=0,
                    help='number of iterations for equilibration')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.02, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=0.2, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=2, help='number of processors')
PARSER.add_argument('--prefix', type=str, default='rpm', help='prefix for all output file names')
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
PARAMS['alpha'] = 6.87098396396261/PARAMS['cubic_side_length']
PARAMS['mu'] = PARAMS['beta_mu']/PARAMS['beta']
PARAMS['charge_plus'] = 1./np.sqrt(1.602176634E-19**2/(4*np.pi*8.8541878128E-12*1e3/1e10/6.02214076E+23))
PARAMS['charge_minus'] = -PARAMS['charge_plus']
PARAMS['dccb_cut'] = PARAMS['cubic_side_length']/int(PARAMS['cubic_side_length']/PARAMS['dccb_cut']) # maximize inside box

def write_feasst_script(params, file_name):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(file_name, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_checkpoint} ln_prob_file {prefix}n{node}_lnpi.txt min_window_size -1
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 1 min_size 3
Checkpoint file_name {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {plus} particle_type1 {minus} cutoff 4.891304347826090 charge0 {charge_plus} charge1 {charge_minus}
Potential VisitModel Ewald alpha {alpha} kmax_squared 38
Potential Model ModelTwoBodyFactory model0 HardSphere model1 ChargeScreened table_size 1e6
ConvertToRefPotential potential_index 1 cutoff {dccb_cut} use_cell true
Potential Model ChargeSelf
ThermoParams beta {beta} chemical_potential0 {mu} chemical_potential1 {mu}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-4

# gcmc initialization and nvt equilibration
TrialAddMultiple particle_type0 0 particle_type1 1 reference_index 0
Log trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_eq.txt
Tune
Run until_num_particles [soft_macro_min] particle_type 0
RemoveTrial name TrialAddMultiple
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} particle_type 0 soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} min_flatness 25 collect_flatness 20 min_collect_sweeps 1
TrialTransferMultiple weight 2 particle_type0 0 particle_type1 1 reference_index 0 num_steps 8
Log trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index].txt
Movie trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_eq.xyz stop_after_iteration 1
Movie trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index].xyz start_after_iteration 1
Tune trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_tune.txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_en.txt multistate true start_after_iteration 1
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_crit.txt
""".format(**params))

def post_process(params):
    lnpi = pd.read_csv(params['prefix']+'n0_lnpi.txt')
    lnpi = lnpi[:3] # cut down to three rows
    lnpi['ln_prob'] -= np.log(sum(np.exp(lnpi['ln_prob'])))  # renormalize
    lnpi['ln_prob_prev'] = [-1.2994315780357, -1.08646312498868, -0.941850889679828]
    lnpi['ln_prob_prev_stdev'] = [0.07, 0.05, 0.05]
    diverged = lnpi[lnpi.ln_prob-lnpi.ln_prob_prev > 5*lnpi.ln_prob_prev_stdev]
    print(diverged)
    assert len(diverged) == 0
    energy = pd.read_csv(params['prefix']+'n0s0_en.txt')
    energy = energy[:3]
    energy['prev'] = [0, -0.939408, -2.02625]
    energy['prev_stdev'] = [1e-14, 0.02, 0.04]
    diverged = energy[energy.average-energy.prev > 10*energy.prev_stdev]
    print(diverged)
    assert len(diverged) == 0

if __name__ == '__main__':
    feasstio.run_simulations(params=PARAMS,
                             sim_node_dependent_params=None,
                             write_feasst_script=write_feasst_script,
                             post_process=post_process,
                             queue_function=feasstio.slurm_single_node,
                             args=ARGS)
