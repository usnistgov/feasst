"""
Tabulate a potential with multiple site types.
At the default temperature, the trimer model, also demonstrated in the flat_histogram tutorial 8,
forms micelles which have difficulty sampling in this canonical ensemble simulation,
especially when there are no cluster moves (see https://doi.org/10.1063/1.4918557).
Regardless, the formation of micelle-like structures demonstrates the tabular potential may be
appropriately accounting for the different site types.
Note that if more than one particle is present, the site types indices add to the ones of the
previous type (e.g., 2 particles with two site types each would have site types 0 and 1 in
particle type 0 and site types 2 and 3 in particle type 1.
"""

import subprocess
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/trimer.txt',
        help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=1./0.2, help='inverse temperature')
    parser.add_argument('--num_particles', type=int, default=500, help='number of particles')
    parser.add_argument('--rho_lower', type=float, default=9e-2, help='lowest number density')
    parser.add_argument('--rho_upper', type=float, default=9e-2, help='highest number density')
    parser.add_argument('--tpc', type=int, default=int(1e5), help='trials per cycle')
    parser.add_argument('--equilibration_cycles', type=int, default=int(1e1),
        help='number of cycles for equilibration')
    parser.add_argument('--production_cycles', type=int, default=int(1e2),
        help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=0.1, help='hours until termination')
    parser.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
    parser.add_argument('--table_file', type=str, default='trimer_table.txt', help='table file name')
    parser.add_argument('--num_z', type=int, default=int(1e3), help='number of table elements')
    parser.add_argument('--inner', type=float, default=0.75, help='As described in TablePotential')
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
    params['prefix'] = 'trimer'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    params['rhos'] = np.linspace(params['rho_lower'], params['rho_upper'], num=params['num_sims'])
    params['cubic_side_lengths'] = np.power(params['num_particles']/params['rhos'], 1./3.).tolist()
    params['rhos'] = params['rhos'].tolist()
    params['gamma'] = -2 # as described in TablePotential
    params['cutoff'] = 3
    params['rwca'] = 2**(1/6.)
    generate_table(params)
    return params, args

def sim_node_dependent_params(params):
    params['cubic_side_length'] = params['cubic_side_lengths'][params['sim']]

def user_potential(distance, site1, site2):
    en_lj = 4*(distance**-12 - distance**-6)
    if site1 == 0 and site2 == 0:
        return en_lj
    else:
        return en_lj + 1

def generate_table(params):
    assert params['num_z'] > 1
    dz = 1./(params['num_z'] - 1)
    with open(params['table_file'], 'w') as file1:
        file1.write("""site_types 2 A R""")
        for site1 in [0, 1]:
            for site2 in [0, 1]:
                if site1 <= site2:
                    if site1 == 0 and site2 == 0:
                        cutoff = params['cutoff']
                    else:
                        cutoff = params['rwca']
                    rcg = cutoff**params['gamma']
                    rhg = params['inner']**params['gamma']
                    file1.write("""\ninner {inner}\nnum_z {num_z}\n""".format(**params))
                    for z in np.arange(0, 1 + dz/2, dz):
                        if z == 0:
                            distance = params['inner']
                        else:
                            distance = (z*(rcg - rhg) + rhg)**(1./params['gamma'])
                        en = user_potential(distance, site1, site2)
                        #print('distance', distance, 'en', en)
                        file1.write(str(en) + " ")

# write fst script to run a simulation
def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed={seed}
Configuration cubic_side_length={cubic_side_length} particle_type=trimer:{fstprt} cutoff={cutoff} cutoffA_R={rwca} cutoffR_R={rwca}
Potential Model=TablePotential table_file={table_file}
Potential VisitModel=LongRangeCorrections
ThermoParams beta={beta} chemical_potential=-1
Metropolis
TrialTranslate tunable_param=2
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# grand canonical ensemble initalization
TrialAdd particle_type=trimer
Run until_num_particles={num_particles}
Remove name=TrialAdd

# canonical ensemble equilibration
Metropolis trials_per_cycle={tpc} cycles_to_complete={equilibration_cycles}
Tune
CheckEnergy trials_per_update={tpc} decimal_places=8
Let [write]=trials_per_write={tpc} output_file={prefix}{sim:03d}
Log [write]_eq.txt
Run until=complete
Remove name=Tune,Log

# canonical ensemble production
Metropolis trials_per_cycle={tpc} cycles_to_complete={production_cycles}
Log [write].csv
Movie [write].xyz
Energy [write]_en.csv
CPUTime [write]_cpu.txt
Run until=complete
""".format(**params))

def post_process(params):
    """ placeholder """
    assert True

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
