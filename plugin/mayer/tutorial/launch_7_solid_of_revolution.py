"""
Example Mayer-sampling simulation of a Spherocylinder.
Compare with https://doi.org/10.1063/1.5004687
For L=1, D=1, B2 = 6.02138591938044
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import accumulator

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--table_file', type=str, default='/feasst/plugin/patch/test/data/tablek5l1.0d1.txt',
                        help='table file name')
    parser.add_argument('--tpc', type=int, default=int(1e4), help='trials per cycle')
    parser.add_argument('--equilibration_cycles', type=int, default=int(1e1),
                        help='number of cycles for equilibration')
    parser.add_argument('--production_cycles', type=int, default=int(1e2),
                        help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=0.1, help='hours until termination')
    parser.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
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
    params['prefix'] = 'solid_rev_'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    return params, args

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length 500 particle_type0 /feasst/plugin/patch/particle/one_patch.txt \
    add_particles_of_type0 2 \
    group0 first first_particle_index 0 \
    group1 centers centers_site_type0 0
Potential Model IdealGas VisitModelInner SolidOfRevolutionTable table_file {table_file} group centers
RefPotential Model HardSphere group centers sigma 0 sigma0 1
ThermoParams beta 1
MayerSampling trials_per_cycle {tpc} cycles_to_complete {equilibration_cycles}
TrialTranslate new_only true reference_index 0 tunable_param 1 group first
TrialRotate new_only true reference_index 0 tunable_param 40
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# tune trial parameters
CriteriaWriter trials_per_write {tpc} output_file {prefix}{sim}_b2_eq.txt
Log trials_per_write {tpc} output_file {prefix}{sim}_eq.txt
Movie trials_per_write {tpc} output_file {prefix}{sim}_eq.xyz
Tune
Run until complete
Remove name0 Tune name1 CriteriaWriter name2 Log name3 Movie

# production
CriteriaWriter trials_per_write {tpc} output_file {prefix}{sim}_b2.txt
Log trials_per_write {tpc} output_file {prefix}{sim}.txt
Movie trials_per_write {tpc} output_file {prefix}{sim}.xyz
MayerSampling trials_per_cycle {tpc} cycles_to_complete {production_cycles}
Run until complete
""".format(**params))

def post_process(params):
    def b2(file_name):
        file1 = open(file_name, 'r')
        lines = file1.readlines()
        file1.close()
        exec('iprm=' + lines[0], globals())
        return iprm
    b2hs_ref = 2*np.pi/3 # reference HS
    b2_overall = accumulator.Accumulator()
    for i in range(params['num_sims']):
        print('i', i)
        b2t = b2(params['prefix']+str(i)+'_b2.txt')
        b2 = b2t['second_virial_ratio']*b2hs_ref
        b2_overall.add(b2)
        print('b2', b2, '+/-', b2t['second_virial_ratio_block_stdev']*b2hs_ref)
        assert np.abs(b2 - 6.02138591938044) < 0.12
    print('b2 overall', b2_overall.mean(), b2_overall.stdev()/np.sqrt(b2_overall.num_values()))

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
