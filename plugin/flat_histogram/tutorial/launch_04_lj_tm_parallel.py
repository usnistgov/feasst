"""
Flat-histogram simulation of single-site Lennard Jones particles in the grand canonical ensemble.
The default temperature is above the critical point.

A common question is "how long does a flat histogram simulation take to finish?"
This tutorially usually finishes in less than an hour with the default parameters.
Because flat histogram is an iterative convergence process, the best way we have to measure its completion is the number of "sweeps" as defined in https://dx.doi.org/10.1063/1.4918557 .

One way to periodically check on the progress of flat histogram simulations is to see the number of sweeps in each window with the following BASH command:
grep num_sweeps lj*_crit.txt

The convergence of the macrostate distrubtion function in the ljn0_lnpi.txt files is also a way to monitor the simulations.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio

def parse(fstprt='/feasst/particle/lj.fstprt',
          beta=1./1.5,
          beta_mu=-1.568214,
          mu_init=10,
          min_sweeps=1e2,
          cubic_side_length=8,
          max_particles=370,
          window_alpha=2.25,
          hours_checkpoint=0.02,
          hours_terminate=0.2,
          num_nodes=1,
          min_window_size=5,
          procs_per_node=32,
          tpc=int(1e6),
          collect_flatness=20,
          min_flatness=25):
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default=fstprt, help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=beta, help='1/(K_Boltzmann*Temperature)')
    parser.add_argument('--beta_mu', type=float, default=beta_mu, help='chemical potential')
    parser.add_argument('--mu_init', type=float, default=mu_init, help='chemical potential to initialize number of particles')
    parser.add_argument('--max_particles', type=int, default=max_particles, help='maximum number of particles')
    parser.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
    parser.add_argument('--window_alpha', type=int, default=window_alpha, help='as window_alpha increases, window size of large N decreases')
    parser.add_argument('--min_window_size', type=int, default=min_window_size, help='minimum window size when parallelizing')
    parser.add_argument('--min_sweeps', type=int, default=min_sweeps,
                        help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
    parser.add_argument('--collect_flatness', type=int, default=collect_flatness, help='WL flatness to begin collection matrix')
    parser.add_argument('--min_flatness', type=int, default=min_flatness, help='WL flatness to begin transition matrix')
    parser.add_argument('--cubic_side_length', type=float, default=cubic_side_length,
                        help='cubic periodic boundary length')
    parser.add_argument('--tpc', type=int, default=tpc, help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=0, help='number of cycles for equilibration')
    parser.add_argument('--hours_checkpoint', type=float, default=hours_checkpoint, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=hours_terminate, help='hours until termination')
    parser.add_argument('--procs_per_node', type=int, default=procs_per_node, help='number of processors')
    parser.add_argument('--run_type', '-r', type=int, default=0,
                        help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--seed', type=int, default=-1,
                        help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
    parser.add_argument('--num_nodes', type=int, default=num_nodes, help='Number of nodes in queue')
    parser.add_argument('--scratch', type=str, default=None,
                        help='Optionally write scheduled job to scratch/logname/jobid.')
    parser.add_argument('--node', type=int, default=0, help='node ID')
    parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
    parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

    # Convert arguments into a parameter dictionary, and add argument-dependent parameters.
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    params = vars(args)
    params['script'] = __file__
    params['prefix'] = 'lj'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.95*params['hours_terminate'] - 0.05 # terminate FEASST before SLURM
    params['hours_terminate'] *= params['procs_per_node'] # real time -> cpu time
    params['hours_checkpoint'] *= params['procs_per_node']
    params['num_sims'] = params['num_nodes']
    params['procs_per_sim'] = params['procs_per_node']
    params['mu'] = params['beta_mu']/params['beta']
    params['system'] = """Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt}
Potential Model LennardJones
Potential VisitModel LongRangeCorrections""".format(**params)
    params['nvt_trials'] = "TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25"
    params['muvt_trials'] = "TrialTransfer weight 2 particle_type 0"
    params['init_trials'] = "TrialAdd particle_type 0"
    params['init_remove'] = "Remove name TrialAdd"
    params['min_particles_second_window'] = ""
    return params, args

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_checkpoint} ln_prob_file {prefix}n{node}_lnpi.txt min_window_size -1
WindowExponential maximum {max_particles} min0 {min_particles} {min_particles_second_window} num {procs_per_node} overlap 0 alpha {window_alpha} min_size {min_window_size}
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

RandomMT19937 seed {seed}
{system}
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
{nvt_trials}
CheckEnergy trials_per_update {tpc} decimal_places 4

# gcmc initialization and nvt equilibration
{init_trials}
Log trials_per_write {tpc} output_file {prefix}n{node}s[sim_index]_eq.txt
Tune
Run until_num_particles [soft_macro_min]
{init_remove}
ThermoParams beta {beta} chemical_potential {mu}
Metropolis trials_per_cycle {tpc} cycles_to_complete {equilibration}
Run until complete
Remove name0 Tune name1 Log

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} min_flatness {min_flatness} collect_flatness {collect_flatness} min_collect_sweeps 1
{muvt_trials}
#To print xyz for each macrostate in separate files, add the following arguments to the "Movie" lines below: multistate true multistate_aggregate false
Movie           trials_per_write {tpc} output_file {prefix}n{node}s[sim_index]_eq.xyz stop_after_cycle 1
Movie           trials_per_write {tpc} output_file {prefix}n{node}s[sim_index].xyz start_after_cycle 1
Log             trials_per_write {tpc} output_file {prefix}n{node}s[sim_index].txt
Tune            trials_per_write {tpc} output_file {prefix}n{node}s[sim_index]_tune.txt multistate true stop_after_cycle 1
Energy          trials_per_write {tpc} output_file {prefix}n{node}s[sim_index]_en.txt multistate true start_after_cycle 1
HeatCapacity    trials_per_write {tpc} output_file {prefix}n{node}s[sim_index]_cv.txt multistate true start_after_cycle 1
CriteriaWriter  trials_per_write {tpc} output_file {prefix}n{node}s[sim_index]_crit.txt
ProfileCPU      trials_per_write {tpc} output_file {prefix}n{node}s[sim_index]_profile.csv
CriteriaUpdater trials_per_update 1e5
""".format(**params))

def post_process(params):
    """ Skip the following checks if temperature is not 1.5 """
    if np.abs(params['beta'] - 1./1.5) > 1e-5:
        return
    lnpi=pd.read_csv(params['prefix']+'n0_lnpi.txt')
    gce_av_num_particles = (np.exp(lnpi["ln_prob"]) * lnpi["state"]).sum()
    assert np.abs(gce_av_num_particles - 310.4179421879679) < 0.5
    srsw = pd.read_csv(params['feasst_install']+'../plugin/flat_histogram/test/data/stat150.csv')
    plt.plot(lnpi['state'], lnpi['ln_prob'], label='FEASST')
    plt.plot(srsw['N'], srsw['lnPI'], linestyle='dashed', label='SRSW')
    plt.xlabel('number of particles', fontsize=16)
    plt.ylabel('ln probability', fontsize=16)
    plt.legend(fontsize=16)
    #plt.savefig(params['prefix']+'_lnpi.png', bbox_inches='tight', transparent='True')

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
