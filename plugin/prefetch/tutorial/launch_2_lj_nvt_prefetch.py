"""
Prefetching canonical ensemble Monte Carlo simulation of Lennard Jones particles.
"""

import argparse
import json
from pyfeasst import fstio
from pyfeasst import physical_constants

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/lj.fstprt',
                        help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=1./0.9, help='inverse temperature')
    parser.add_argument('--num_particles', type=int, default=400, help='number of particles')
    parser.add_argument('--cubic_side_length', type=float, default=8,
                        help='cubic periodic boundary length')
    parser.add_argument('--tpc', type=int, default=int(1e4), help='trials per cycle')
    parser.add_argument('--equilibration_cycles', type=int, default=int(1e2),
                        help='number of cycles for equilibraiton')
    parser.add_argument('--production_cycles', type=int, default=int(1e3),
                        help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
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
    params['prefix'] = 'lj'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_node'] = 1
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    params['target_acceptance'] = 0.25
    if params['fstprt'] == '/feasst/particle/lj.fstprt':
        params['potential'] = """Potential Model LennardJones
    Potential VisitModel LongRangeCorrections"""
        params['trials'] = """TrialTranslate tunable_param 2 tunable_target_acceptance {target_acceptance}""".format(**params)
    elif params['fstprt'] == '/feasst/particle/spce.fstprt':
        params['cubic_side_length'] = 20
        params['num_particles'] = 265
        params['beta'] = 1./(300*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
        params['alpha'] = 5.6/params['cubic_side_length']
        params['dccb_cut'] = 0.9*3.165
        params['dccb_cut'] = params['cubic_side_length']/int(params['cubic_side_length']/params['dccb_cut']) # maximize inside box
        params['potential'] = """Potential VisitModel Ewald alpha {alpha} kmax_squared 38
    Potential Model ModelTwoBodyFactory model0 LennardJones model1 ChargeScreened erfc_table_size 2e4 VisitModel VisitModelCutoffOuter
    #RefPotential Model HardSphere group oxygen cutoff {dccb_cut} VisitModel VisitModelCell min_length {dccb_cut} cell_group oxygen
    Potential Model ChargeScreenedIntra VisitModel VisitModelBond
    Potential Model ChargeSelf
    Potential VisitModel LongRangeCorrections""".format(**params)
        params['trials'] = """TrialTranslate tunable_param 0.2 tunable_target_acceptance {target_acceptance}
    TrialParticlePivot weight 0.5 particle_type 0 tunable_param 0.5 tunable_target_acceptance {target_acceptance} """.format(**params)
    else:
        assert False, 'unrecognized fstprt'
    return params, args

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
CheckEnergy trials_per_update {tpc} decimal_places 4

# gcmc initialization
TrialAdd particle_type 0
Log trials_per_write {tpc} output_file {prefix}{sim}_init.txt
Tune
Run until_num_particles {num_particles}
Remove name0 TrialAdd name1 Log

# nvt equilibration
ThermoParams beta {beta}
Metropolis trials_per_cycle {tpc} cycles_to_complete {equilibration_cycles}
Log trials_per_write {tpc} output_file {prefix}{sim}_eq.txt
Movie trials_per_write {tpc} output_file {prefix}{sim}_eq.xyz
Run until complete
Remove name0 Tune name1 Log name2 Movie

# nvt production
Metropolis trials_per_cycle {tpc} cycles_to_complete {production_cycles}
Log trials_per_write {tpc} output_file {prefix}{sim}.txt
Movie trials_per_write {tpc} output_file {prefix}{sim}.xyz
Energy trials_per_write {tpc} output_file {prefix}{sim}_en.txt append true
CPUTime trials_per_write {tpc} output_file {prefix}{sim}_time.txt append true
Run until complete
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
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
