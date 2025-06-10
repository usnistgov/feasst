"""
Post process a simulation by stepping through the configurations of an xyz file.
Note that the xyz file given has few configurations to save space in the git repo,
and thus the results are very noisey.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import scattering

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/lj_new.txt',
                        help='FEASST particle definition')
    parser.add_argument('--xyz_file', type=str, default='post_process.xyz',
                        help='The xyz file to read configurations')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=1., help='hours until termination')
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
    params['prefix'] = 'post_process'
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
Configuration particle_type=fluid:{fstprt} xyz_file={xyz_file} cutoff=1
Potential Model=HardSphere VisitModel=VisitModelCell min_length=max_cutoff
ThermoParams beta=1 chemical_potential=1
Metropolis
Scattering trials_per_update=1 trials_per_write=1 num_frequency=10 output_file={prefix}{sim:03d}_iq.csv
ReadConfigFromFile input_file={xyz_file}
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}
Run until=complete
""".format(**params))

def post_process(params):
    iq=pd.read_csv(params['prefix']+'000_iq.csv', comment="#")
    grp = iq.groupby('q', as_index=False)
    assert np.abs(iq['i'][3810] - 0.114066) < 0.4
    assert np.abs(iq['i'][0]/iq['p0'][0]**2 - 1) < 0.075
    plt.scatter(iq['q'], iq['i']/iq['p0']**2, label='sq_all')
    plt.plot(grp.mean()['q'], grp.mean()['i']/grp.mean()['p0']**2, label='sq_av')
    plt.legend()
    #plt.savefig(params['prefix']+'.png', bbox_inches='tight', transparent='True')

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
