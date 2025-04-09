"""
Flat-histogram simulation of monovalent RPM in the grand canonical expanded ensemble.
Using expanded ensemble, single particle insertions and deletions are possible.
Dual-cut configurational bias is used for insertions and deletions.
The results are are compared with https://doi.org/10.1063/1.5123683
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import physical_constants
from pyfeasst import macrostate_distribution

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--plus', type=str, default='/feasst/plugin/charge/particle/rpm_plus.fstprt',
                        help='FEASST particle definition of the positive charge RPM')
    parser.add_argument('--minus', type=str, default='/feasst/plugin/charge/particle/rpm_minus.fstprt',
                        help='FEASST particle definition of the negative charge RPM')
    parser.add_argument('--beta', type=float, default=1./0.047899460618081, help='inverse temperature')
    parser.add_argument('--beta_mu', type=float, default=-13.94, help='beta times chemical potential')
    parser.add_argument('--max_particles', type=int, default=10, help='maximum number of particles')
    parser.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
    parser.add_argument('--min_sweeps', type=int, default=1e2,
                        help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
    parser.add_argument('--cubic_side_length', type=float, default=12,
                        help='cubic periodic boundary length')
    parser.add_argument('--dccb_cut', type=float, default=2**(1./6.),
                        help='dual-cut configurational bias cutoff')
    parser.add_argument('--tpc', type=int, default=int(1e6), help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=0,
                        help='number of cycles for equilibration')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
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
    params['prefix'] = 'rpm'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.95*params['hours_terminate'] - 0.05 # terminate FEASST before SLURM
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    params['alpha'] = 6.87098396396261/params['cubic_side_length']
    params['mu'] = params['beta_mu']/params['beta']
    params['charge_plus'] = 1./np.sqrt(1.602176634E-19**2/(4*np.pi*8.8541878128E-12*1e3/1e10/6.02214076E+23))
    params['charge_minus'] = -params['charge_plus']
    params['dccb_cut'] = params['cubic_side_length']/int(params['cubic_side_length']/params['dccb_cut']) # maximize inside box
    params['minimums'] = [params['min_particles']]
    if params['num_nodes'] == 1:
        params['windows'] = macrostate_distribution.window_exponential(
            alpha=2.5, minimums=params['minimums'], maximum=params['max_particles'],
            number=params['num_sims'], overlap=1, min_size=2)
    return params, args

def sim_node_dependent_params(params):
    """ Define parameters that are dependent on the sim or node. """
    params['min_particles'] = params['windows'][params['sim']][0]
    params['max_particles'] = params['windows'][params['sim']][1]
    params['sim_start'] = 0
    params['sim_end'] = params['num_sims'] - 1

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {plus} particle_type1 {minus} cutoff 4.891304347826090 charge0 {charge_plus} charge1 {charge_minus}
Potential VisitModel Ewald alpha {alpha} kmax_squared 38
Potential Model ModelTwoBodyFactory model0 HardSphere model1 ChargeScreened erfc_table_size 2e4
RefPotential Model HardSphere VisitModel VisitModelCell min_length {dccb_cut}
Potential Model ChargeSelf
ThermoParams beta {beta} chemical_potential0 {mu} chemical_potential1 {mu}
Metropolis Constraint AEqualB extra_A 1
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
CheckEnergy trials_per_update {tpc} decimal_places 4
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# gcmc initialization and nvt equilibration
TrialAdd weight 1 particle_type 0 reference_index 0 num_steps 8
TrialAdd weight 1 particle_type 1 reference_index 0 num_steps 8
Log trials_per_write {tpc} output_file {prefix}n{node}s{sim}_eq.csv
Tune
Run until_num_particles {min_particles}
Remove name0 TrialAdd name1 TrialAdd
Metropolis trials_per_cycle {tpc} cycles_to_complete {equilibration}
Run until complete
Remove name0 Tune name1 Log

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} \
  Bias WLTM min_sweeps {min_sweeps} min_flatness 25 collect_flatness 20 min_collect_sweeps 1 \
  Constraint AEqualB extra_A 1
TrialTransfer weight 1 particle_type 0 reference_index 0 num_steps 8
TrialTransfer weight 1 particle_type 1 reference_index 0 num_steps 8
Log            trials_per_write {tpc} output_file {prefix}n{node}s{sim}.csv
Movie          trials_per_write {tpc} output_file {prefix}n{node}s{sim}_eq.xyz stop_after_cycle 1
Movie          trials_per_write {tpc} output_file {prefix}n{node}s{sim}.xyz start_after_cycle 1
Tune           trials_per_write {tpc} output_file {prefix}n{node}s{sim}_tune.csv multistate true stop_after_cycle 1
Energy         trials_per_write {tpc} output_file {prefix}n{node}s{sim}_en.csv multistate true start_after_cycle 1
CriteriaWriter trials_per_write {tpc} output_file {prefix}n{node}s{sim:03d}_crit.csv
CriteriaUpdater trials_per_update 1e5
Run until complete

# continue until all simulations on the node are complete
WriteFileAndCheck sim {sim} sim_start {sim_start} sim_end {sim_end} file_prefix {prefix}n{node}s file_suffix _finished.txt output_file {prefix}n{node}_terminate.txt
Run until_file_exists {prefix}n{node}_terminate.txt trials_per_file_check {tpc}
""".format(**params))

def post_process(params):
    lnpi=macrostate_distribution.splice_files(prefix=params['prefix']+'n0s', suffix='_crit.csv', shift=False).dataframe()
    lnpi=lnpi[:5] # cut down to five rows
    lnpi=lnpi[::2] # drop odd rows for neutral-only macrostates
    lnpi['ln_prob'] -= np.log(sum(np.exp(lnpi['ln_prob'])))  # renormalize
    lnpi['ln_prob_prev'] = [-1.2994315780357, -1.08646312498868, -0.941850889679828]
    lnpi['ln_prob_prev_stdev'] = [0.07, 0.05, 0.05]
    diverged = lnpi[lnpi.ln_prob-lnpi.ln_prob_prev > 5*lnpi.ln_prob_prev_stdev]
    print(diverged)
    assert len(diverged) == 0
    energy = pd.read_csv(params['prefix']+'n0s0_en.csv')
    energy = energy[:5]
    energy['prev'] = [0, -0.115474, -0.939408, -1.32485, -2.02625]
    energy['prev_stdev'] = [1e-14, 1e-6, 0.02, 0.03, 0.04]
    diverged = energy[energy.average-energy.prev > 10*energy.prev_stdev]
    print(diverged)
    assert len(diverged) == 0

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
