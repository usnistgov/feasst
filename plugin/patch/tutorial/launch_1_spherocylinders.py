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
    parser.add_argument('--num_particles', type=int, default=500, help='number of particles')
    parser.add_argument('--cylinder_length', type=float, default=3, help='cylinder length (distance between center of end caps)')
    parser.add_argument('--cubic_side_length', type=int, default=20, help='cubic periodic boundary length')
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
    params['prefix'] = 'spherocylinders'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_node'] = 1
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']

    hard_spherocylinder(length=params['cylinder_length'],
                        diameter=1,
                        file_name=params['prefix']+'.txt')
    return params, args

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed={seed}
Configuration cubic_side_length={cubic_side_length} particle_type=cylinder:{prefix}.txt group=centers centers_site_type=0
Potential Model=HardSphere VisitModelInner=Spherocylinder group=centers
ThermoParams beta=1 chemical_potential=-1
Metropolis
TrialTranslate tunable_param=2
TrialRotate tunable_param=40
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}
CheckEnergy trials_per_update={tpc} decimal_places=6

# gcmc initialization
TrialAdd particle_type=cylinder
Let [write]=trials_per_write={tpc} output_file={prefix}{sim:03d}
Log [write]_init.csv
Tune
Run until_num_particles={num_particles}
Remove name=TrialAdd,Log

# nvt equilibration
Metropolis trials_per_cycle={tpc} cycles_to_complete={equilibration}
Log [write]_eq.csv
Movie [write]_eq.xyz
Run until=complete
Remove name=Tune,Log,Movie

# nvt production
Metropolis trials_per_cycle={tpc} cycles_to_complete={production}
Log [write].csv
Movie [write].xyz
MovieSpherocylinder [write]c.xyz
Energy [write]_en.csv
Run until=complete
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
