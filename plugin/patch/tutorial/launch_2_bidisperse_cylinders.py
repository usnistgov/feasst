"""
Simulate hard spherocylinders.
Distances are based on the diameter of the cylinder (e.g., unit diameter).
"""

import argparse
import json
from pyfeasst import fstio
from make_spherocylinder import hard_spherocylinder

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    #parser.add_argument('--fstprt', type=str, default='/feasst/plugin/patch/particle/spherocylinder.txt',
    #                    help='FEASST particle definition')
    parser.add_argument('--num_particles', type=int, default=100, help='number of particles')
    parser.add_argument('--small_cylinder_length', type=float, default=2, help='small cylinder length (distance between center of end caps)')
    parser.add_argument('--small_cylinder_diameter', type=float, default=1, help='small cylinder diameter')
    parser.add_argument('--large_cylinder_length', type=float, default=10, help='large cylinder length (distance between center of end caps)')
    parser.add_argument('--large_cylinder_diameter', type=float, default=3, help='large cylinder diameter')
    parser.add_argument('--cubic_side_length', type=int, default=30, help='cubic periodic boundary length')
    parser.add_argument('--tpc', type=int, default=int(1e4), help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=int(1e1),
                        help='number of cycles for equilibraiton')
    parser.add_argument('--production', type=int, default=int(1e2),
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
    params['prefix'] = 'bidisp_cyl'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_node'] = 1
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']

    for size in ['small', 'large']:
        hard_spherocylinder(length=params[size+'_cylinder_length'],
                            diameter=params[size+'_cylinder_diameter'],
                            file_name=params['prefix'] + '_' + size + '.txt')

    # write initial config of two parallel large cylinders
    with open(params['prefix'] + '_init.xyz', 'w') as f:
        xyz_params = {'cubic_side_length': params['cubic_side_length'],
                      'half_separation': 2}
        f.write("""4
-1 {cubic_side_length} {cubic_side_length} {cubic_side_length} 0 0 0
0 -{half_separation} 0 0
1 -{half_separation} 0 1
0  {half_separation} 0 0
1  {half_separation} 0 1""".format(**xyz_params))
    return params, args

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {prefix}_small.txt \
    particle_type1 {prefix}_large.txt \
    group0 centers centers_site_type0 0 centers_site_type1 2 \
    add_particles_of_type1 2 xyz_file {prefix}_init.xyz
Potential Model HardSphere VisitModelInner Spherocylinder group centers
ThermoParams beta 1 chemical_potential0 -1 chemical_potential1 -1
Metropolis
TrialTranslate weight 1 tunable_param 2 particle_type 0
TrialTranslate weight 1 tunable_param 0.1 particle_type 1 dimension 0
TrialRotate weight 1 tunable_param 40 particle_type 0
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}
CheckEnergy trials_per_update {tpc} tolerance 1e-4

# gcmc initialization
TrialAdd particle_type 0
Log trials_per_write {tpc} output_file {prefix}{sim}_init.txt
Tune
Run until_num_particles {num_particles}
Remove name0 TrialAdd name1 Log

# nvt equilibration
Metropolis trials_per_cycle {tpc} cycles_to_complete {equilibration}
Log trials_per_write {tpc} output_file {prefix}{sim}_eq.txt
Movie trials_per_write {tpc} output_file {prefix}{sim}_eq.xyz
Run until complete
Remove name0 Tune name1 Log name2 Movie

# nvt production
Metropolis trials_per_cycle {tpc} cycles_to_complete {production}
Log trials_per_write {tpc} output_file {prefix}{sim}.txt
Movie trials_per_write {tpc} output_file {prefix}{sim}.xyz
MovieSpherocylinder trials_per_write {tpc} output_file {prefix}{sim}c.xyz
Energy trials_per_write {tpc} output_file {prefix}{sim}_en.txt
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
