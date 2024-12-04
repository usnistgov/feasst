"""
Simulate a mixture of CO2 and N2.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import macrostate_distribution
from pyfeasst import physical_constants

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt0', type=str, default='/feasst/particle/co2.fstprt',
                    help='FEASST particle definition of the first particle.')
PARSER.add_argument('--fstprt1', type=str, default='/feasst/particle/n2.fstprt',
                    help='FEASST particle definition of the second particle.')
PARSER.add_argument('--temperature', type=float, default=300, help='temperature in Kelvin')
PARSER.add_argument('--mu0', type=float, default=-15.24, help='chemical potential')
PARSER.add_argument('--mu1', type=float, default=-15.24, help='chemical potential')
PARSER.add_argument('--mu_init', type=float, default=10, help='initial chemical potential')
PARSER.add_argument('--num_particles', type=int, default=10, help='total number of particles')
PARSER.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
PARSER.add_argument('--min_sweeps', type=int, default=1e2,
                    help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
PARSER.add_argument('--cubic_side_length', type=float, default=30,
                    help='cubic periodic boundary length')
PARSER.add_argument('--tpc', type=int, default=int(1e6), help='trials per cycle')
PARSER.add_argument('--equilibration_cycles', type=int, default=0,
                    help='number of cycles for equilibration')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.02, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=0.2, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=2, help='number of processors')
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
PARAMS['prefix'] = 'co2n2_'
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.95*PARAMS['hours_terminate'] - 0.05 # terminate FEASST before SLURM
PARAMS['hours_terminate'] *= PARAMS['procs_per_node'] # real time -> cpu time
PARAMS['hours_checkpoint'] *= PARAMS['procs_per_node']
PARAMS['num_sims'] = PARAMS['num_nodes']
PARAMS['procs_per_sim'] = PARAMS['procs_per_node']
PARAMS['ewald_alpha'] = 5.6/PARAMS['cubic_side_length']
PARAMS['beta'] = 1./(PARAMS['temperature']*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_checkpoint} ln_prob_file {prefix}n{node}_lnpi.txt min_window_size -1
WindowExponential maximum {num_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 1.1 min_size 5
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt0} particle_type1 {fstprt1}
Potential VisitModel Ewald alpha {ewald_alpha} kmax_squared 38
Potential Model ModelTwoBodyFactory model0 LennardJones model1 ChargeScreened erfc_table_size 2e4 VisitModel VisitModelCutoffOuter
Potential Model ChargeScreenedIntra VisitModel VisitModelBond
Potential Model ChargeSelf
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential0 {mu_init} chemical_potential1 {mu_init}
Metropolis
TrialTranslate weight 0.5 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 0.5 particle_type 0 tunable_param 0.5 tunable_target_acceptance 0.25
CheckEnergy trials_per_update {tpc} tolerance 1e-4

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Log trials_per_write {tpc} output_file {prefix}n{node}s[sim_index]_eq.txt
Tune
Run until_num_particles [soft_macro_min]
Remove name TrialAdd
TrialAdd particle_type 1
Run until_num_particles {num_particles}
Remove name TrialAdd
ThermoParams beta {beta} chemical_potential0 {mu0} chemical_potential1 {mu1}
Metropolis trials_per_cycle {tpc} cycles_to_complete {equilibration_cycles}
Run until complete
Remove name0 Tune name1 Log

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles particle_type 0 width 1 max {num_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} min_flatness 25 collect_flatness 20 min_collect_sweeps 1
TrialMorph particle_type0 0 particle_type_morph0 1
TrialMorph particle_type0 1 particle_type_morph0 0
Log trials_per_write {tpc} output_file {prefix}n{node}s[sim_index].txt
Movie trials_per_write {tpc} output_file {prefix}n{node}s[sim_index]_eq.xyz stop_after_cycle 1
Movie trials_per_write {tpc} output_file {prefix}n{node}s[sim_index].xyz start_after_cycle 1
Tune trials_per_write {tpc} output_file {prefix}n{node}s[sim_index]_tune.txt multistate true stop_after_cycle 1
Energy trials_per_write {tpc} output_file {prefix}n{node}s[sim_index]_en.txt multistate true start_after_cycle 1
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {tpc} output_file {prefix}n{node}s[sim_index]_crit.txt
""".format(**params))

def post_process(params):
    lnpi = macrostate_distribution.MacrostateDistribution(file_name=params['prefix']+'n0_lnpi.txt')
    #lnpi.plot(show=True)
    assert np.abs(5.137717334901432 - lnpi.average_macrostate()) < 0.5

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
