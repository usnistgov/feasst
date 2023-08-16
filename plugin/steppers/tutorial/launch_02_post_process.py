"""
Post process a simulation by stepping through the configurations of an xyz file.
Note that the xyz file given has few configurations to save space in the git repo,
and thus the results are very noisey.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import feasstio
from pyfeasst import scattering

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/forcefield/atom.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--xyz_file', type=str, default='post_process.xyz',
                    help='The xyz file to read configurations')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.2, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=1., help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
PARSER.add_argument('--prefix', type=str, default='post_process', help='prefix for all output file names')
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
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.99*PARAMS['hours_terminate'] - 0.0333 # terminate before queue
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']

def write_feasst_script(params, file_name):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(file_name, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
Configuration particle_type0 {fstprt} xyz_file {xyz_file}
Potential Model HardSphere VisitModel VisitModelCell min_length 1
ThermoParams beta 1 chemical_potential 1
Metropolis
Scattering trials_per_update 1 trials_per_write 1 num_frequency 10 file_name {prefix}{sim}_iq.csv
ReadConfigFromFile file_name {xyz_file}
Checkpoint file_name {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}
Run until_criteria_complete true
""".format(**params))

def post_process(params):
    iq=pd.read_csv(params['prefix']+'0_iq.csv', comment="#")
    grp = iq.groupby('q', as_index=False)
    assert np.abs(iq['i'][3810] - 6.5772) < 0.4
    assert np.abs(iq['i'][0]/iq['p0'][0]**2 - 1) < 0.075
    plt.scatter(iq['q'], iq['i']/iq['p0']**2, label='sq_all')
    plt.plot(grp.mean()['q'], grp.mean()['i']/grp.mean()['p0']**2, label='sq_av')
    plt.legend()
    #plt.savefig(params['prefix']+'.png', bbox_inches='tight', transparent='True')

if __name__ == '__main__':
    feasstio.run_simulations(params=PARAMS,
                             sim_node_dependent_params=None,
                             write_feasst_script=write_feasst_script,
                             post_process=post_process,
                             queue_function=feasstio.slurm_single_node,
                             args=ARGS)
