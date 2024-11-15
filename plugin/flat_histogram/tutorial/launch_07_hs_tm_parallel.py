"""
Flat-histogram simulation of single-site hard sphere particles in the grand canonical ensemble.
Compare equation of state in SRSW:
https://www.nist.gov/mml/csd/chemical-informatics-research-group/hard-sphere-thermodynamic-and-transport-properties
"""

import argparse
import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import macrostate_distribution

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/atom.fstprt',
                        help='FEASST particle definition')
    parser.add_argument('--mu', type=float, default=-2.352321, help='chemical potential')
    parser.add_argument('--mu_init', type=float, default=10, help='initial chemical potential')
    parser.add_argument('--max_particles', type=int, default=256, help='maximum number of particles')
    parser.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
    parser.add_argument('--min_sweeps', type=int, default=1e2,
                        help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
    parser.add_argument('--cubic_side_length', type=float, default=8,
                        help='cubic periodic boundary length')
    parser.add_argument('--tpi', type=int, default=int(1e6),
                        help='trials per iteration, similar to MC cycles, but not necessary num_particles')
    parser.add_argument('--equilibration_iterations', type=int, default=0,
                        help='number of iterations for equilibration')
    parser.add_argument('--hours_checkpoint', type=float, default=0.02, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=0.2, help='hours until termination')
    parser.add_argument('--procs_per_node', type=int, default=32, help='number of processors')
    parser.add_argument('--run_type', '-r', type=int, default=0,
                        help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--seed', type=int, default=-1,
                        help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
    parser.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
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
    params['prefix'] = 'hs'
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
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 1.5 min_size 5
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt}
Potential Model HardSphere VisitModel VisitModelCell min_length 1
ThermoParams beta 1 chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
CheckEnergy trials_per_update {tpi} decimal_places 4

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Log trials_per_write {tpi} output_file {prefix}n{node}s[sim_index]_eq.txt
Tune
Run until_num_particles [soft_macro_min]
RemoveTrial name TrialAdd
ThermoParams beta 1 chemical_potential {mu}
Metropolis num_trials_per_iteration {tpi} num_iterations_to_complete {equilibration_iterations}
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} min_flatness 25 collect_flatness 20 min_collect_sweeps 1
TrialTransfer weight 2 particle_type 0
Log              trials_per_write {tpi} output_file {prefix}n{node}s[sim_index].txt
Movie            trials_per_write {tpi} output_file {prefix}n{node}s[sim_index]_eq.xyz stop_after_iteration 1
Movie            trials_per_write {tpi} output_file {prefix}n{node}s[sim_index].xyz start_after_iteration 1
Tune             trials_per_write {tpi} output_file {prefix}n{node}s[sim_index]_tune.txt multistate true stop_after_iteration 1
Energy           trials_per_write {tpi} output_file {prefix}n{node}s[sim_index]_en.txt multistate true start_after_iteration 1
CriteriaWriter   trials_per_write {tpi} output_file {prefix}n{node}s[sim_index]_crit.txt
PairDistribution trials_per_write {tpi} output_file {prefix}n{node}s[sim_index]_gr.txt \
  dr 0.025 multistate true multistate_aggregate false trials_per_update 1000
CriteriaUpdater trials_per_update 1e5
""".format(**params))

def post_process(params):
    # compare to EOS in SRSW: https://www.nist.gov/mml/csd/chemical-informatics-research-group/hard-sphere-thermodynamic-and-transport-properties
    lnpi = macrostate_distribution.MacrostateDistribution(file_name=params['prefix']+'n0_lnpi.txt')
    volume = 8**3
    srsw = pd.read_csv(params['feasst_install']+'../plugin/flat_histogram/test/data/stat_hs.csv')
    srsw = srsw[:6]
    pressure = list()
    lnpi_rw = copy.deepcopy(lnpi)
    for target_density in srsw['dens']:
        lnpi_rw.reweight_to_macrostate(target_macrostate=target_density*volume)
        pressure.append(-lnpi_rw.ln_prob().values[0]/volume)
    srsw['P_FST'] = pressure
    print(srsw[['dens', 'P_MC', 'P_FST', '+/-']])
    assert np.any(abs(srsw['P_MC'] - srsw['P_FST']) < 1e-3)

    # Use chemical potential from Carnahan-Starling to compare expected average density
    # http://www.sklogwiki.org/SklogWiki/index.php/Carnahan-Starling_equation_of_state
    rho = 0.1
    cubic_side_length = 8
    eta = np.pi/6*rho
    betamu_ex = (8*eta-9*eta**2+3*eta**3)/(1-eta)**3
    betamu = betamu_ex + np.log(rho)
    lnpi_rw = lnpi.reweight(delta_beta_mu=betamu+2.352321)
    density = lnpi_rw.average_macrostate()/volume
    print('target_density', rho, 'density', density)
    assert abs(rho - density) < 5e-4

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
