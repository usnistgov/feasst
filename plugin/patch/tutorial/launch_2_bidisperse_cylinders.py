"""
Simulate hard spherocylinders.
Distances are based on the diameter of the cylinder (e.g., unit diameter).
"""

import argparse
import json
from pyfeasst import fstio
from make_spherocylinder import hard_spherocylinder

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
#PARSER.add_argument('--fstprt', type=str, default='/feasst/plugin/patch/particle/spherocylinder.fstprt',
#                    help='FEASST particle definition')
PARSER.add_argument('--num_particles', type=int, default=100, help='number of particles')
PARSER.add_argument('--small_cylinder_length', type=float, default=2, help='small cylinder length (distance between center of end caps)')
PARSER.add_argument('--small_cylinder_diameter', type=float, default=1, help='small cylinder diameter')
PARSER.add_argument('--large_cylinder_length', type=float, default=10, help='large cylinder length (distance between center of end caps)')
PARSER.add_argument('--large_cylinder_diameter', type=float, default=3, help='large cylinder diameter')
PARSER.add_argument('--cubic_side_length', type=int, default=30, help='cubic periodic boundary length')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e4),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=int(1e1),
                    help='number of iterations for equilibraiton')
PARSER.add_argument('--production_iterations', type=int, default=int(1e2),
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
PARAMS['prefix'] = 'bidisp_cyl'
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.99*PARAMS['hours_terminate'] - 0.0333 # terminate before queue
PARAMS['procs_per_node'] = 1
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']

for size in ['small', 'large']:
    hard_spherocylinder(length=PARAMS[size+'_cylinder_length'],
                        diameter=PARAMS[size+'_cylinder_diameter'],
                        file_name=PARAMS['prefix'] + '_' + size + '.fstprt')

# write initial config of two parallel large cylinders
with open(PARAMS['prefix'] + '_init.xyz', 'w') as f:
    xyz_params = {'cubic_side_length': PARAMS['cubic_side_length'],
                  'half_separation': 2}
    f.write("""4
-1 {cubic_side_length} {cubic_side_length} {cubic_side_length} 0 0 0
0 -{half_separation} 0 0
1 -{half_separation} 0 1
0  {half_separation} 0 0
1  {half_separation} 0 1""".format(**xyz_params))

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {prefix}_small.fstprt \
    particle_type1 {prefix}_large.fstprt \
    group0 centers centers_site_type0 0 centers_site_type1 2 \
    add_particles_of_type1 2 xyz_file {prefix}_init.xyz
Potential Model HardSphere VisitModelInner Spherocylinder group centers
ThermoParams beta 1 chemical_potential0 -1 chemical_potential1 -1
Metropolis
TrialTranslate weight 1 tunable_param 2 particle_type 0
TrialTranslate weight 1 tunable_param 0.1 particle_type 1 dimension 0
TrialRotate weight 1 tunable_param 40 particle_type 0
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-4

# gcmc initialization
TrialAdd particle_type 0
Log trials_per_write {trials_per_iteration} output_file {prefix}{sim}_init.txt
Tune
Run until_num_particles {num_particles}
RemoveTrial name TrialAdd
RemoveAnalyze name Log

# nvt equilibration
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Log trials_per_write {trials_per_iteration} output_file {prefix}{sim}_eq.txt
Movie trials_per_write {trials_per_iteration} output_file {prefix}{sim}_eq.xyz
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log
RemoveAnalyze name Movie

# nvt production
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {production_iterations}
Log trials_per_write {trials_per_iteration} output_file {prefix}{sim}.txt
Movie trials_per_write {trials_per_iteration} output_file {prefix}{sim}.xyz start_after_iteration 1
MovieSpherocylinder trials_per_write {trials_per_iteration} output_file {prefix}{sim}c.xyz start_after_iteration 1
Energy trials_per_write {trials_per_iteration} output_file {prefix}{sim}_en.txt
Run until_criteria_complete true
""".format(**params))

def post_process(params):
    assert True

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
