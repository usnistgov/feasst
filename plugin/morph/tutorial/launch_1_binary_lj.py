"""
Simulate Mixture I. of https://doi.org/10.1063/1.1844372
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import macrostate_distribution

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt0', type=str, default='/feasst/particle/atom.fstprt',
                        help='FEASST particle definition of the first particle.')
    parser.add_argument('--fstprt1', type=str, default='/feasst/particle/lj.fstprt',
                        help='FEASST particle definition of the second particle.')
    parser.add_argument('--beta', type=float, default=0.8, help='inverse temperature')
    parser.add_argument('--mu0', type=float, default=-5.4, help='chemical potential')
    parser.add_argument('--mu1', type=float, default=-5.5, help='chemical potential')
    parser.add_argument('--mu_init', type=float, default=10, help='initial chemical potential')
    parser.add_argument('--num_particles', type=int, default=20, help='total number of particles')
    parser.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
    parser.add_argument('--min_sweeps', type=int, default=1e2,
                        help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
    parser.add_argument('--cubic_side_length', type=float, default=7,
                        help='cubic periodic boundary length')
    parser.add_argument('--tpc', type=int, default=int(1e6), help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=0,
                        help='number of cycles for equilibration')
    parser.add_argument('--hours_checkpoint', type=float, default=0.02, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=0.2, help='hours until termination')
    parser.add_argument('--procs_per_node', type=int, default=2, help='number of processors')
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
    params['prefix'] = 'binary_lj'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.95*params['hours_terminate'] - 0.05 # terminate FEASST before SLURM
    params['hours_terminate'] *= params['procs_per_node'] # real time -> cpu time
    params['hours_checkpoint'] *= params['procs_per_node']
    params['num_sims'] = params['num_nodes']
    params['procs_per_sim'] = params['procs_per_node']
    return params, args

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_checkpoint} ln_prob_file {prefix}n{node}_lnpi.txt min_window_size -1
WindowExponential maximum {num_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 1.1 min_size 5
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt0} particle_type1 {fstprt1} sigma0 1.0 epsilon0 1.0 cutoff0 3.0 sigma1 1.064 epsilon1 1.37 cutoff1 3.0 sigma0_1 1.034 epsilon0_1 1.152 cutoff0_1 3.0
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential0 {mu_init} chemical_potential1 {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
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
Metropolis trials_per_cycle {tpc} cycles_to_complete {equilibration}
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
    assert np.abs(9.327631384282558 - lnpi.average_macrostate()) < 0.5

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
