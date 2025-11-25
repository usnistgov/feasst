"""
This is a simple test of a ModelTableSphere1D by computing the Henry
Coefficient of a particle in a sphere.
"""

import argparse
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/lj_new.txt',
                        help='FEASST particle definition')
    parser.add_argument('--plot_table', type=int, default=0, help='0: no plot, 1: plot')
    parser.add_argument('--radius', type=float, default=3, help='radius of cylindrical cavity')
    parser.add_argument('--beta', type=float, default=1., help='inverse temperature')
    parser.add_argument('--mu', type=float, default=-1, help='chemical potential')
    parser.add_argument('--cubic_side_length', type=float, default=9,
                        help='cubic periodic boundary length')
    parser.add_argument('--tpc', type=int, default=int(1e5), help='trials per cycle')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
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
    params['prefix'] = 'sph1d'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.95*params['hours_terminate'] - 0.05 # terminate FEASST before SLURM
    params['num_sims'] = params['num_nodes']
    params['procs_per_sim'] = params['procs_per_node']
    with open(params['prefix']+'_shape.txt', 'w') as file1:
        file1.write("Sphere radius=3 center=0,0,0")
    def user_potential(r, # distance from inner surface, not radial distance from center
                       radius):
        return 1./np.power(r, 3) - 1./np.power(radius, 3)
    with open(params['prefix']+'_table.txt', 'w') as file1:
        file1.write("site_types=LJ\n")
        radius = 3
        dr = 0.01
        r = np.arange(dr, radius+dr/2, 0.01)
        file1.write("999999999 ")
        values = user_potential(r, radius)
        for v in values:
            file1.write(str(v)+" ")
        file1.write("\n")
    if args.plot_table == 1:
        # check the energy interpolated from the table against the analytical value
        dists = np.arange(0.5, params['radius'], 0.033)
        ens = list()
        for dist in dists:
            params['dist'] = 3-dist
            run_en(params)
            df = pd.read_csv(params['prefix']+'_one.csv')
            #ens.append(df['ModelLJShape'].values[0])
            ens.append(df['ModelTableSphere1D'].values[0])
        import matplotlib.pyplot as plt
        plt.plot(dists, ens, label='table')
        plt.plot(dists, user_potential(dists, radius), color='black', linestyle='dotted', label='analytical')
        plt.xlabel('r', fontsize=16)
        plt.ylabel('U', fontsize=16)
        plt.legend(fontsize=16)
        plt.show()
        quit()
    return params, args

def run_en(params):
    """
    Run a simulation to obtain the energy of a particle in confinement using params['dist'].
    """
    with open(params['prefix']+"_one.xyz", "w") as file1: file1.write(
"""1
-1 8 8 8
0 {dist} 0 0""".format(**params))
    with open(params['prefix']+"_launch.txt", "w") as myfile: myfile.write("""
MonteCarlo
RandomMT19937 seed={seed}
Configuration xyz_file={prefix}_one.xyz particle_type={fstprt}
#Potential Model=ModelLJShape shape_file={prefix}_shape.txt
Potential Model=ModelTableSphere1D radius=3 center=0,0,0 table_file={prefix}_table.txt
ThermoParams beta=1000000
Metropolis
Log output_file={prefix}_one.csv max_precision=true clear_file=true
Run num_trials=1
""".format(**params))
    syscode = subprocess.call(params['feasst_install']+"""bin/fst < {prefix}_launch.txt > {prefix}_launch.log""".format(**params), shell=True, executable='/bin/bash')
    if syscode > 0: sys.exit(1)

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed={seed}
Configuration cubic_side_length={cubic_side_length} particle_type=fluid:{fstprt}
Potential Model=HardSphere
#Potential Model=ModelLJShape shape_file={prefix}_shape.txt
Potential Model=ModelTableSphere1D radius=3 center=0,0,0 table_file={prefix}_table.txt
ThermoParams beta={beta} chemical_potential={mu}
AlwaysReject
TrialAdd particle_type=fluid new_only=true
HenryCoefficient trials_per_write={tpc} output_file={prefix}.csv write_precision=12 num_beta_taylor=4
Run num_trials=1e6
""".format(**params))

def post_process(params):
    df = pd.read_csv(params['prefix'] + '.csv', comment='#')
    print(df)
    assert np.abs(df['average'][0] - 0.03764) < 3*df['block_stdev'][0]
    with open(params['prefix']+'.csv') as f:
        firstline = f.readline().rstrip()
        henry=eval(firstline[1:])
        print(henry)
        assert np.abs(henry['beta_taylor'][0] - 0.03764) < 3*df['block_stdev'][0]
        assert np.abs(henry['beta_taylor'][1]) < 0.1

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
