"""
Canonical ensemble Monte Carlo simulation of Lennard Jones particles in a triclinic domain with cell lists.
VisitModelCell in triclinic domains increases min_length by the ratio of shorted length to largest inscribed sphere.
If min_length is too small to ensure all interactions are in neighboring cells, then the CheckEnergy will produce an error if OptPotential is used and not equivalent to the given Potential.
"""

import argparse
import json
from pyfeasst import fstio

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/lj.fstprt',
                        help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=1./1.5, help='inverse temperature')
    parser.add_argument('--num_particles', type=int, default=500, help='number of particles')
    parser.add_argument('--cubic_side_length', type=float, default=20, help='cubic periodic boundary length')
    parser.add_argument('--tpc', type=int, default=int(1e4), help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=int(1e1), help='number of cycles for equilibraiton')
    parser.add_argument('--production', type=int, default=int(1e2), help='number of cycles for production')
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
    params['prefix'] = 'triclinic'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_node'] = 1
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    return params, args

def sim_node_dependent_params(params):
    """ Set parameters that depent upon the sim or node here. """

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt} xy 5 xz 5 yz 5
Potential Model LennardJones
OptPotential Model LennardJones VisitModel VisitModelCell min_length 3
Potential VisitModel LongRangeCorrections
OptPotential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential -1
Metropolis
TrialTranslate weight 1 tunable_param 2 tunable_target_acceptance 2
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}
CheckEnergy trials_per_update {tpc} decimal_places 4

# gcmc initialization
TrialAdd particle_type 0
Log trials_per_write {tpc} output_file {prefix}{sim}_init.csv
Tune
Run until_num_particles {num_particles}
Remove name0 TrialAdd name1 Log name2 Tune

# nvt equilibration
ThermoParams beta {beta}
Metropolis trials_per_cycle {tpc} cycles_to_complete {equilibration}
Tune trials_per_tune 20
Log     trials_per_write {tpc} output_file {prefix}{sim}_eq.csv
Movie   trials_per_write {tpc} output_file {prefix}{sim}_eq.xyz
Run until complete
Remove name0 Tune name1 Log name2 Movie

# nvt production
Metropolis trials_per_cycle {tpc} cycles_to_complete {production}
Log        trials_per_write {tpc} output_file {prefix}{sim}.csv
Movie      trials_per_write {tpc} output_file {prefix}{sim}.xyz
Energy     trials_per_write {tpc} output_file {prefix}{sim}_en.csv
ProfileCPU trials_per_write {tpc} output_file {prefix}{sim}_profile.csv
Run until complete
""".format(**params))

def post_process(params):
    """ place holder """

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
