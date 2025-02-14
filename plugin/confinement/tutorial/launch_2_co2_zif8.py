"""
Flat-histogram simulation of TraPPE CO2 adsorption in ZIF8.
Compare with https://doi.org/10.1021/jp400480q .
The next tutorial with bulk CO2 is required to obtain the pressure.
The ZIF8 forcefield is described in https://doi.org/10.1002/chem.200902144 .
Developed with Dr. Siderius in 2024 and maintained by Dr. Hatch.
"""

import argparse
import subprocess
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
    parser.add_argument('--fluid', type=str, default='/feasst/particle/co2.fstprt',
                        help='FEASST particle definition of a fluid particle / adsorbate.')
    parser.add_argument('--MOF', type=str, default='/feasst/plugin/confinement/particle/ZIF8_rep222_PerezPellitero.fstprt',
                        help='FEASST particle definition of the MOF / adsorbent.')
    parser.add_argument('--cutoff', type=float, default=15, help='site-site cutoff distance in Angstroms')
    parser.add_argument('--temperature', type=float, default=303, help='temperature in Kelvin')
    parser.add_argument('--mu', type=float, default=-20, help='chemical potential')
    parser.add_argument('--mu_init', type=float, default=-1, help='initial chemical potential')
    parser.add_argument('--beta_init', type=float, default=0.1, help='initial beta to fill adsorbent')
    parser.add_argument('--max_particles', type=int, default=150, help='maximum number of particles')
    parser.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
    parser.add_argument('--min_sweeps', type=int, default=1,
                        help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
    parser.add_argument('--cubic_side_length', type=float, default=34.023240,
                        help='cubic periodic boundary length')
    parser.add_argument('--tpc', type=int, default=int(1e5), help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=1e0, help='number of cycles for equilibration')
    parser.add_argument('--hours_checkpoint', type=float, default=0.5, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
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
    params['prefix'] = 'co2-zif8_'
    params['alpha'] = 5.6/params['cubic_side_length']
    params['beta'] = 1./(params['temperature']*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.95*params['hours_terminate'] - 0.05 # terminate FEASST before SLURM
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    params['min1'] = params['min_particles'] + 10
    params['minimums'] = [params['min_particles'], params['min1']]
    if params['num_nodes'] == 1:
        params['windows'] = macrostate_distribution.window_exponential(
            alpha=2.2, minimums=params['minimums'], maximum=params['max_particles'],
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
Configuration cubic_side_length {cubic_side_length} particle_type0 {fluid} particle_type1 {MOF} add_particles_of_type1 1 group0 fluid fluid_particle_type 0 group1 MOF MOF_particle_type 1 cutoff {cutoff}
Potential VisitModel Ewald alpha {alpha} kmax_squared 27
Potential Model ModelTwoBodyFactory model0 LennardJonesForceShift model1 ChargeScreened erfc_table_size 2e4 VisitModel VisitModelCutoffOuter
Potential Model ChargeScreenedIntra VisitModel VisitModelBond
Potential Model ChargeSelf
ZeroBackground
ThermoParams beta {beta_init} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 30 particle_type 0 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 20 particle_type 0 tunable_param 0.5 tunable_target_acceptance 0.25
CheckEnergy trials_per_update {tpc} tolerance 1e-4
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# write the MOF xyz for visualization
Movie output_file {prefix}n{node}s{sim}_MOF.xyz group MOF clear_file true
WriteStepper analyze_name Movie
Remove name Movie

# gcmc initialization and nvt equilibration
TrialAdd weight 50.0 particle_type 0
Log trials_per_write {tpc} output_file {prefix}n{node}s{sim}_eq.csv
Tune
Run until_num_particles {min_particles} particle_type 0
Remove name TrialAdd
ThermoParams beta {beta} chemical_potential {mu}
Metropolis trials_per_cycle {tpc} cycles_to_complete {equilibration}
Run until complete
Remove name0 Tune name1 Log

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles particle_type 0 width 1 max {max_particles} min {min_particles} \
  Bias TransitionMatrix min_sweeps {min_sweeps}
TrialTransfer weight 50 particle_type 0
Log trials_per_write {tpc} output_file {prefix}n{node}s{sim}.csv
#To print trajectories for each macrostate in separate files, add the following arguments to the "Movie" lines below: multistate true multistate_aggregate false
Movie trials_per_write {tpc} output_file {prefix}n{node}s{sim}_eq.xyz stop_after_cycle 1 group fluid
Movie trials_per_write {tpc} output_file {prefix}n{node}s{sim}.xyz start_after_cycle 1 group fluid
Tune trials_per_write {tpc} output_file {prefix}n{node}s{sim}_tune.csv multistate true stop_after_cycle 1
Energy trials_per_write {tpc} output_file {prefix}n{node}s{sim}_en.csv multistate true start_after_cycle 1
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {tpc} output_file {prefix}n{node}s{sim}_crit.csv
Run until complete

# continue until all simulations on the node are complete
WriteFileAndCheck sim {sim} sim_start {sim_start} sim_end {sim_end} file_prefix {prefix}n{node}s file_suffix _finished.txt output_file {prefix}n{node}_terminate.txt
Run until_file_exists {prefix}n{node}_terminate.txt trials_per_file_check {tpc}
""".format(**params))

def post_process(params):
    assert True # placeholder
    print('launching after_2_co2_bulk.py')
    subprocess.check_call("""python after_2_co2_bulk.py --run_type {run_type}""".format(**params), shell=True, executable='/bin/bash')

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
