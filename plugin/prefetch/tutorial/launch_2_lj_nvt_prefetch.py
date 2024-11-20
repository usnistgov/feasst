"""
Prefetching canonical ensemble Monte Carlo simulation of Lennard Jones particles.
"""

import argparse
import json
from pyfeasst import fstio
from pyfeasst import physical_constants

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/particle/lj.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--beta', type=float, default=1./0.9, help='inverse temperature')
PARSER.add_argument('--num_particles', type=int, default=400, help='number of particles')
PARSER.add_argument('--cubic_side_length', type=float, default=8,
                    help='cubic periodic boundary length')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e4),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=int(1e2),
                    help='number of iterations for equilibraiton')
PARSER.add_argument('--production_iterations', type=int, default=int(1e3),
                    help='number of iterations for production')
PARSER.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
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
PARAMS['procs_per_node'] = 1
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']
PARAMS['target_acceptance'] = 0.25
if PARAMS['fstprt'] == '/feasst/particle/lj.fstprt':
    PARAMS['potential'] = """Potential Model LennardJones
Potential VisitModel LongRangeCorrections"""
    PARAMS['trials'] = """TrialTranslate tunable_param 2 tunable_target_acceptance {target_acceptance}""".format(**PARAMS)
elif PARAMS['fstprt'] == '/feasst/particle/spce.fstprt':
    PARAMS['cubic_side_length'] = 20
    PARAMS['num_particles'] = 265
    PARAMS['beta'] = 1./(300*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
    PARAMS['alpha'] = 5.6/PARAMS['cubic_side_length']
    PARAMS['dccb_cut'] = 0.9*3.165
    PARAMS['dccb_cut'] = PARAMS['cubic_side_length']/int(PARAMS['cubic_side_length']/PARAMS['dccb_cut']) # maximize inside box
    PARAMS['potential'] = """Potential VisitModel Ewald alpha {alpha} kmax_squared 38
Potential Model ModelTwoBodyFactory model0 LennardJones model1 ChargeScreened erfc_table_size 2e4 VisitModel VisitModelCutoffOuter
#RefPotential Model HardSphere group oxygen cutoff {dccb_cut} VisitModel VisitModelCell min_length {dccb_cut} cell_group oxygen
Potential Model ChargeScreenedIntra VisitModel VisitModelBond
Potential Model ChargeSelf
Potential VisitModel LongRangeCorrections""".format(**PARAMS)
    PARAMS['trials'] = """TrialTranslate tunable_param 0.2 tunable_target_acceptance {target_acceptance}
TrialParticlePivot weight 0.5 particle_type 0 tunable_param 0.5 tunable_target_acceptance {target_acceptance} """.format(**PARAMS)
else:
    assert False, 'unrecognized fstprt'

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
Prefetch synchronize true
#Prefetch synchronize true trials_per_check 1
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt}
{potential}
ThermoParams beta {beta} chemical_potential -1
Metropolis
{trials}
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}
CheckEnergy trials_per_update {trials_per_iteration} decimal_places 4

# gcmc initialization
TrialAdd particle_type 0
Log trials_per_write {trials_per_iteration} output_file {prefix}{sim}_init.txt
Tune
Run until_num_particles {num_particles}
Remove name0 TrialAdd name1 Log

# nvt equilibration
ThermoParams beta {beta}
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Log trials_per_write {trials_per_iteration} output_file {prefix}{sim}_eq.txt
Movie trials_per_write {trials_per_iteration} output_file {prefix}{sim}_eq.xyz
Run until_criteria_complete true
Remove name0 Tune name1 Log name2 Movie

# nvt production
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {production_iterations}
Log trials_per_write {trials_per_iteration} output_file {prefix}{sim}.txt
Movie trials_per_write {trials_per_iteration} output_file {prefix}{sim}.xyz start_after_iteration 1
Energy trials_per_write {trials_per_iteration} output_file {prefix}{sim}_en.txt append true
CPUTime trials_per_write {trials_per_iteration} output_file {prefix}{sim}_time.txt append true
Run until_criteria_complete true
""".format(**params))

def linear_fit(x, b):
    return -0.5*x + b

def post_process(params):
    """ Compute standard deviation of energy with time """
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.optimize import curve_fit
    time = pd.read_csv('lj0_time.txt', delim_whitespace=True, header=None)
    print(time[1])
    en = pd.read_csv('lj0_en.txt', header=None, comment="a")
    print(en[0])
    equil=500
    logt = np.log(time[1][equil:])
    logs = np.log(en[2][equil:])
    popt, pcov = curve_fit(linear_fit, logt, logs)
    plt.scatter(np.log(time[1]), np.log(en[2]))
    plt.plot(logt, linear_fit(logt, popt[0]), color='black')
    print(popt[0])
    #plt.xscale('log')
    #plt.yscale('log')
    plt.show()

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
