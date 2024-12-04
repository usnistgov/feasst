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
PARSER.add_argument('--num_particles', type=int, default=500, help='number of particles')
PARSER.add_argument('--cylinder_length', type=float, default=3, help='cylinder length (distance between center of end caps)')
PARSER.add_argument('--cubic_side_length', type=int, default=20, help='cubic periodic boundary length')
PARSER.add_argument('--tpc', type=int, default=int(1e4), help='trials per cycle')
PARSER.add_argument('--equilibration', type=int, default=int(1e1),
                    help='number of cycles for equilibraiton')
PARSER.add_argument('--production', type=int, default=int(1e2),
                    help='number of cycles for production')
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
PARAMS['prefix'] = 'spherocylinders'
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.99*PARAMS['hours_terminate'] - 0.0333 # terminate before queue
PARAMS['procs_per_node'] = 1
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']

hard_spherocylinder(length=PARAMS['cylinder_length'],
                    diameter=1,
                    file_name=PARAMS['prefix']+'.fstprt')

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {prefix}.fstprt group0 centers centers_site_type0 0
Potential Model HardSphere VisitModelInner Spherocylinder group centers
ThermoParams beta 1 chemical_potential -1
Metropolis
TrialTranslate weight 1 tunable_param 2
TrialRotate weight 1 tunable_param 40
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
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
