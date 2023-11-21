"""
Example single-site Lennard-Jones Gibbs ensemble Monte Carlo simulation using FEASST.
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
PARSER.add_argument('--beta', type=float, default=1./0.8, help='inverse temperature')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e5),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=int(1e1),
                    help='number of iterations for equilibraiton')
PARSER.add_argument('--production_iterations', type=int, default=int(5e2),
                    help='number of iterations for production')
PARSER.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
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
PARAMS['prefix'] = 'lj'
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.99*PARAMS['hours_terminate'] - 0.0333 # terminate before queue
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']
def sim_node_dependent_params(params):
    """ Set parameters that depent upon the sim or node here. """

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length 12 particle_type0 {fstprt}
Configuration cubic_side_length 8 particle_type0 {fstprt}
Potential Model LennardJones configuration_index 0
Potential Model LennardJones configuration_index 1
#Potential Model LennardJones VisitModel VisitModelCell min_length max_cutoff configuration_index 0
#Potential Model LennardJones VisitModel VisitModelCell min_length max_cutoff configuration_index 1
Potential VisitModel LongRangeCorrections configuration_index 0
Potential VisitModel LongRangeCorrections configuration_index 1
RefPotential VisitModel DontVisitModel configuration_index 0
RefPotential VisitModel DontVisitModel configuration_index 1
ThermoParams beta {beta} chemical_potential 10
Metropolis
TrialTranslate tunable_param 2 tunable_target_acceptance 0.2 configuration_index 0
TrialTranslate tunable_param 0.1 tunable_target_acceptance 0.2 configuration_index 1
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# grand canonical ensemble initalization
Log trials_per_write {trials_per_iteration} output_file {prefix}{sim}_fill.csv
Movie trials_per_write {trials_per_iteration} output_file {prefix}{sim}_c0_fill.xyz configuration_index 0
Movie trials_per_write {trials_per_iteration} output_file {prefix}{sim}_c1_fill.xyz configuration_index 1
Tune
TrialAdd particle_type 0 configuration_index 0
Run until_num_particles 32 configuration_index 0
RemoveTrial name TrialAdd
TrialAdd particle_type 0 configuration_index 1
Run until_num_particles 480 configuration_index 1
RemoveTrial name TrialAdd
RemoveAnalyze name Log
RemoveAnalyze name Movie
RemoveAnalyze name Movie

# gibbs ensemble equilibration
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
TrialGibbsParticleTransfer weight 0.05 particle_type 0 reference_index 0
TrialGibbsVolumeTransfer weight 0.001 tunable_param 0.1 reference_index 0
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-8
CheckConstantVolume trials_per_update {trials_per_iteration} tolerance 1e-4
Log trials_per_write {trials_per_iteration} output_file {prefix}{sim}_eq.csv
Movie trials_per_write {trials_per_iteration} output_file {prefix}{sim}_c0_eq.xyz configuration_index 0
Movie trials_per_write {trials_per_iteration} output_file {prefix}{sim}_c1_eq.xyz configuration_index 1
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log
RemoveAnalyze name Movie
RemoveAnalyze name Movie

# gibbs ensemble production
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {production_iterations}
Log trials_per_write {trials_per_iteration} output_file {prefix}{sim}.csv
Movie trials_per_write {trials_per_iteration} output_file {prefix}{sim}_c0.xyz configuration_index 0
Movie trials_per_write {trials_per_iteration} output_file {prefix}{sim}_c1.xyz configuration_index 1
Energy trials_per_write {trials_per_iteration} output_file {prefix}{sim}_c0_en.csv configuration_index 0
Energy trials_per_write {trials_per_iteration} output_file {prefix}{sim}_c1_en.csv configuration_index 1
Density trials_per_write {trials_per_iteration} output_file {prefix}{sim}_c0_dens.csv configuration_index 0
Density trials_per_write {trials_per_iteration} output_file {prefix}{sim}_c1_dens.csv configuration_index 1
PressureFromTestVolume trials_per_update 1e3 trials_per_write {trials_per_iteration} output_file {prefix}{sim}_pressure.csv
# the pressure in the liquid phase is harder to converge with test volume changes? and faster to compute in the vapor.
#PressureFromTestVolume trials_per_update 1e3 trials_per_write {trials_per_iteration} output_file {prefix}{sim}_c1_pressure.csv configuration_index 1
CPUTime trials_per_write {trials_per_iteration} output_file {prefix}{sim}_cpu.csv
Run until_criteria_complete true
""".format(**params))

def post_process(params):
    z_factor = 3
    vapor_density = pd.read_csv(params['prefix']+"0_c0_dens.csv")
    vapor_density['diff'] = np.abs(vapor_density['average']-6.1007E-03)
    vapor_density['tol'] = np.sqrt(vapor_density['block_stdev']**2+5.63E-07**2)
    #print(vapor_density)
    diverged = vapor_density[vapor_density['diff'] > z_factor*vapor_density['tol']]
    print(diverged)
    #assert len(diverged) == 0
    liquid_density = pd.read_csv(params['prefix']+"0_c1_dens.csv")
    liquid_density['diff'] = np.abs(liquid_density['average']-0.79981)
    liquid_density['tol'] = np.sqrt(liquid_density['block_stdev']**2+0.000013**2)
    #print(liquid_density)
    diverged = liquid_density[liquid_density['diff'] > z_factor*liquid_density['tol']]
    print(diverged)
    #assert len(diverged) == 0
    pressure = pd.read_csv(params['prefix']+"0_pressure.csv")
    pressure['diff'] = np.abs(pressure['pressure_average']-0.0046465)
    pressure['tol'] = np.sqrt(pressure['pressure_block_stdev']**2+(3.74e-7)**2)
    #print(pressure)
    diverged = pressure[pressure['diff'] > z_factor*pressure['tol']]
    print(diverged)
    #assert len(diverged) == 0

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
