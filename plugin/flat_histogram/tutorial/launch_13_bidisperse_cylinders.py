"""
Simulate two large hard spherocylinders in a bath of small hard spherocylinders.
Use flat histogram to compute the depletant effect as a function of separation distance
along the first dimension direction.
"""

import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
sys.path.insert(0, '../../patch/tutorial/')
from make_spherocylinder import hard_spherocylinder

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/lj.fstprt',
                        help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=1./1.5, help='inverse temperature')
    parser.add_argument('--num_particles', type=int, default=100, help='number of particles')
    parser.add_argument('--small_cylinder_length', type=float, default=0, help='small cylinder length (distance between center of end caps)')
    parser.add_argument('--small_cylinder_diameter', type=float, default=1, help='small cylinder diameter')
    parser.add_argument('--large_cylinder_length', type=float, default=0, help='large cylinder length (distance between center of end caps)')
    parser.add_argument('--large_cylinder_diameter', type=float, default=3, help='large cylinder diameter')
    parser.add_argument('--displacement', type=float, default=0.1, help='displacement steps of large cylinder')
    parser.add_argument('--mu_init', type=float, default=10, help='initial chemical potential')
    parser.add_argument('--min_sweeps', type=int, default=1,
                        help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
    parser.add_argument('--cubic_side_length', type=float, default=10,
                        help='cubic periodic boundary length')
    parser.add_argument('--tpc', type=int, default=int(1e6), help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=0, help='number of cycles for equilibration')
    parser.add_argument('--hours_checkpoint', type=float, default=0.02, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=0.2, help='hours until termination')
    parser.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
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
    params['prefix'] = 'bidisp_sphc'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.95*params['hours_terminate'] - 0.05 # terminate FEASST before SLURM
    params['hours_terminate'] *= params['procs_per_node'] # real time -> cpu time
    params['hours_checkpoint'] *= params['procs_per_node']
    params['num_sims'] = params['num_nodes']
    params['procs_per_sim'] = params['procs_per_node']

    # write initial config of two parallel large cylinders
    with open(params['prefix'] + '_init.xyz', 'w') as f:
        xyz_params = {'cubic_side_length': params['cubic_side_length'],
                      'half_separation': 2}
        f.write("""4
    -1 {cubic_side_length} {cubic_side_length} {cubic_side_length} 0 0 0
    0 {half_separation} 0 0
    1 {half_separation} 0 1
    0 -{half_separation} 0 0
    1 -{half_separation} 0 1""".format(**xyz_params))

    for size in ['small', 'large']:
        hard_spherocylinder(length=params[size+'_cylinder_length'],
                            diameter=params[size+'_cylinder_diameter'],
                            file_name=params['prefix'] + '_' + size + '.fstprt')
    return params, args

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {prefix}_small.fstprt \
    particle_type1 {prefix}_large.fstprt \
    particle_type2 {prefix}_large.fstprt \
    group0 centers centers_site_type0 0 centers_site_type1 2 centers_site_type2 4 \
    add_particles_of_type1 1 add_particles_of_type2 1 xyz_file {prefix}_init.xyz
Potential Model HardSphere VisitModelInner Spherocylinder group centers
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 2 particle_type 0
TrialRotate weight 1 tunable_param 40 particle_type 0
CheckEnergy trials_per_update {tpc} decimal_places 4

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Log trials_per_write {tpc} output_file {prefix}n{node}_eq.txt
Tune
Run until_num_particles {num_particles}
Remove name TrialAdd
Metropolis trials_per_cycle {tpc} cycles_to_complete {equilibration}
Run until complete
Remove name0 Tune name1 Log

# gcmc tm production
FlatHistogram Macrostate MacrostatePosition particle_index 0 site_index 0 dimension 0 width 0.1 max 2.1001 min 1.1001 \
  Bias WLTM min_sweeps {min_sweeps} min_flatness 25 collect_flatness 20 min_collect_sweeps 1
TrialTranslate weight 1 tunable_param 0.1 particle_type 1 dimension 0
Log                 trials_per_write {tpc} output_file {prefix}n{node}.txt
Movie               trials_per_write {tpc} output_file {prefix}n{node}_eq.xyz stop_after_cycle 1
MovieSpherocylinder trials_per_write {tpc} output_file {prefix}n{node}_eqc.xyz stop_after_cycle 1
Movie               trials_per_write {tpc} output_file {prefix}n{node}.xyz start_after_cycle 1
MovieSpherocylinder trials_per_write {tpc} output_file {prefix}n{node}c.xyz start_after_cycle 1
Tune                trials_per_write {tpc} output_file {prefix}n{node}_tune.txt multistate true stop_after_cycle 1
Energy              trials_per_write {tpc} output_file {prefix}n{node}_en.txt multistate true start_after_cycle 1
CriteriaWriter      trials_per_write {tpc} output_file {prefix}n{node}_crit.txt
CriteriaUpdater     trials_per_update 1e5
Run until complete
""".format(**params))

def post_process(params):
    assert True

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
