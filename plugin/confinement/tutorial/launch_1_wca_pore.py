"""
Flat-histogram simulation of single-site Lennard Jones particles in the grand canonical ensemble.
Simulate the adsorption of an LJ fluid instead of a rigid porous network of WCA particles,
as described in https://doi.org/10.1021/acs.jpcb.3c00613 .
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
PARSER.add_argument('--pore', type=str, default='/feasst/plugin/confinement/particle/porel4.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--beta', type=float, default=1./0.8, help='inverse temperature')
PARSER.add_argument('--mu', type=float, default=-1, help='chemical potential')
PARSER.add_argument('--mu_init', type=float, default=10, help='initial chemical potential')
PARSER.add_argument('--max_particles', type=int, default=64, help='maximum number of particles')
PARSER.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
PARSER.add_argument('--min_sweeps', type=int, default=2,
                    help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
PARSER.add_argument('--cubic_side_length', type=float, default=9,
                    help='cubic periodic boundary length')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e5),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=1e0,
                    help='number of iterations for equilibration')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.02, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=0.2, help='hours until termination')
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
PARAMS['prefix'] = 'pore'
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.95*PARAMS['hours_terminate'] - 0.05 # terminate FEASST before SLURM
PARAMS['hours_terminate'] *= PARAMS['procs_per_node'] # real time -> cpu time
PARAMS['hours_checkpoint'] *= PARAMS['procs_per_node']
PARAMS['num_sims'] = PARAMS['num_nodes']
PARAMS['procs_per_sim'] = PARAMS['procs_per_node']
PARAMS['wca'] = 2**(1./6.)
PARAMS['dccb_cut'] = PARAMS['cubic_side_length']/int(PARAMS['cubic_side_length']/PARAMS['wca']) # maximize inside box
PARAMS['min1'] = PARAMS['min_particles'] + 2

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_checkpoint} ln_prob_file {prefix}{node}_lnpi.txt min_window_size -1
WindowExponential maximum {max_particles} min0 {min_particles} min1 {min1} num {procs_per_node} overlap 0 alpha 1.0 min_size 2
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt} particle_type1 {pore} cutoff0 2.5 cutoff1 {wca} cutoff0_1 {wca} add_particles_of_type1 1 group0 liquid liquid_particle_type 0 group1 pore pore_particle_type 1
NeighborCriteria maximum_distance 1.375 minimum_distance 0.9 site_type0 0 site_type1 0
Potential EnergyMap EnergyMapNeighborCriteria neighbor_index 0 Model LennardJonesForceShift
ConvertToRefPotential cutoff {dccb_cut} use_cell true
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
TrialAVB2 weight 0.1 particle_type 0
TrialAVB4 weight 0.1 particle_type 0
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-4

# write the pore xyz for visualization
Movie output_file {prefix}{node}s[sim_index]_pore.xyz group pore clear_file true
WriteStepper analyze_name Movie
RemoveAnalyze name Movie

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Log trials_per_write {trials_per_iteration} output_file {prefix}{node}s[sim_index]_eq.txt
Tune
Run until_num_particles [soft_macro_min] particle_type 0
RemoveTrial name TrialAdd
ThermoParams beta {beta} chemical_potential {mu}
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles particle_type 0 width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} min_flatness 25 collect_flatness 20 min_collect_sweeps 1
TrialTransfer weight 2 particle_type 0 reference_index 0 num_steps 4
TrialTransferAVB weight 0.2 particle_type 0 reference_index 0 num_steps 4
Log trials_per_write {trials_per_iteration} output_file {prefix}{node}s[sim_index].txt
#To print trajectories for each macrostate in separate files, add the following arguments to the "Movie" lines below: multistate true multistate_aggregate false
Movie trials_per_write {trials_per_iteration} output_file {prefix}{node}s[sim_index]_eq.xyz stop_after_iteration 1 group liquid
Movie trials_per_write {trials_per_iteration} output_file {prefix}{node}s[sim_index].xyz start_after_iteration 1 group liquid
Tune trials_per_write {trials_per_iteration} output_file {prefix}{node}s[sim_index]_tune.txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per_iteration} output_file {prefix}{node}s[sim_index]_en.txt multistate true start_after_iteration 1
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {trials_per_iteration} output_file {prefix}{node}s[sim_index]_crit.txt
""".format(**params))

def post_process(params):
    import pandas as pd
    df = pd.read_csv(params['prefix']+'0_lnpi.txt')
    #print(df['delta_ln_prob'])
    delta_ln_prob_1_2_3 = [5.33710, 4.66512, 4.27458]
    #print(df['delta_ln_prob_stdev'])
    delta_ln_prob_std_1_2_3 = [0.002927, 0.006865, 0.005270]
    df = df[1:4]
    df['delta_ln_prob_prev'] = [5.33710, 4.66512, 4.27458]
    df['delta_ln_prob_std_prev'] = [0.002927, 0.006865, 0.005270]
    print(df)
    z_factor = 6
    diff = df[df.delta_ln_prob - df.delta_ln_prob_prev > z_factor*df.delta_ln_prob_stdev]
    print(diff)
    assert len(diff) == 0

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
