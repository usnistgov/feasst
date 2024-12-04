"""
This tutorial is essentially a copy of /path/to/feasst/tutorial/launch.py,
except that the LJ potential is tabulated.
This table potential is nearly as fast as the built-in LJ potential,
and this tutorial serves a test of the table potential.
"""

import subprocess
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/particle/lj.fstprt',
    help='FEASST particle definition')
PARSER.add_argument('--beta', type=float, default=1./0.9, help='inverse temperature')
PARSER.add_argument('--num_particles', type=int, default=500, help='number of particles')
PARSER.add_argument('--rho_lower', type=float, default=1e-3, help='lowest number density')
PARSER.add_argument('--rho_upper', type=float, default=9e-3, help='highest number density')
PARSER.add_argument('--tpc', type=int, default=int(1e5), help='trials per cycle')
PARSER.add_argument('--equilibration_cycles', type=int, default=int(1e1),
    help='number of cycles for equilibration')
PARSER.add_argument('--production_cycles', type=int, default=int(1e3),
    help='number of cycles for production')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.1, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=0.1, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=5, help='number of processors')
PARSER.add_argument('--table_file', type=str, default='lj_table.txt', help='table file name')
PARSER.add_argument('--plot_table', type=int, default=0, help='0: no plot, 1: plot')
PARSER.add_argument('--num_z', type=int, default=int(1e3), help='number of table elements')
PARSER.add_argument('--inner', type=float, default=0.75, help='As described in TablePotential')
PARSER.add_argument('--cutoff', type=float, default=3, help='potential cutoff distance')
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
PARAMS['prefix'] = 'lj'
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.99*PARAMS['hours_terminate'] - 0.0333 # terminate before queue
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']
PARAMS['rhos'] = np.linspace(PARAMS['rho_lower'], PARAMS['rho_upper'], num=PARAMS['num_sims'])
PARAMS['cubic_side_lengths'] = np.power(PARAMS['num_particles']/PARAMS['rhos'], 1./3.).tolist()
PARAMS['rhos'] = PARAMS['rhos'].tolist()
PARAMS['gamma'] = -2 # as described in TablePotential
def sim_node_dependent_params(params):
    params['cubic_side_length'] = params['cubic_side_lengths'][params['sim']]

def user_potential(distance):
    return 4*(distance**-12 - distance**-6)

def generate_table(PARAMS):
    assert PARAMS['num_z'] > 1
    dz = 1./(PARAMS['num_z'] - 1)
    rhg = PARAMS['inner']**PARAMS['gamma']
    rcg = PARAMS['cutoff']**PARAMS['gamma']
    with open(PARAMS['table_file'], 'w') as file1:
        file1.write("""site_types 1 0\ninner {inner}\nnum_z {num_z}\n""".format(**PARAMS))
        #file1.write("""site_types 1 0\ngamma {gamma}\ninner {inner}\nnum_z {num_z}\n""".format(**PARAMS))
        for z in np.arange(0, 1 + dz/2, dz):
            if z == 0:
                distance = PARAMS['inner']
            else:
                distance = (z*(rcg - rhg) + rhg)**(1./PARAMS['gamma'])
            en = user_potential(distance)
            #print('distance', distance, 'en', en)
            file1.write(str(en) + " ")
generate_table(PARAMS)

def run_en():
    """
    Run a feasst simulation to obtain the energy between two particles as a
    function of sepration distance (params['displacement_test'])
    """
    with open("lj_two.xyz", "w") as file1: file1.write(
"""2
-1 8 8 8
0 0 0 0
1 0 0 {displacement_test}""".format(**PARAMS))
    with open("launch.txt", "w") as myfile: myfile.write("""
MonteCarlo
RandomMT19937 seed time
Configuration xyz_file lj_two.xyz particle_type0 {fstprt} cutoff {cutoff}
Potential Model TablePotential table_file {table_file}
Potential VisitModel LongRangeCorrections
ThermoParams beta 1000000
Metropolis
Log output_file lj.csv max_precision true clear_file true
Run num_trials 1
""".format(**PARAMS))
    syscode = subprocess.call(PARAMS['feasst_install']+"bin/fst < launch.txt > launch.log", shell=True, executable='/bin/bash')
    if syscode > 0: sys.exit(1)

if ARGS.plot_table == 1:
    # check the energy interpolated from the table against the analytical value
    dists = np.arange(0.97, PARAMS['cutoff'], 0.01)
    ens = list()
    for dist in dists:
        PARAMS['displacement_test'] = dist
        run_en()
        df = pd.read_csv('lj.csv')
        ens.append(df['TablePotential'].values[0])
    import matplotlib.pyplot as plt
    plt.plot(dists, ens, label='table')
    plt.plot(dists, user_potential(dists), color='black', linestyle='dotted', label='analytical')
    plt.xlabel('r', fontsize=16)
    plt.ylabel('U', fontsize=16)
    plt.legend(fontsize=16)
    plt.show()

# write fst script to run a simulation
def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt}
Potential Model TablePotential table_file {table_file} VisitModel VisitModelCell min_length max_cutoff
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential -1
Metropolis
TrialTranslate tunable_param 2 tunable_target_acceptance 0.2
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# grand canonical ensemble initalization
TrialAdd particle_type 0
Run until_num_particles {num_particles}
Remove name TrialAdd

# canonical ensemble equilibration
Metropolis trials_per_cycle {tpc} cycles_to_complete {equilibration_cycles}
Tune
CheckEnergy trials_per_update {tpc} tolerance 1e-8
Log trials_per_write {tpc} output_file {prefix}{sim}_eq.txt
Run until complete
Remove name0 Tune name1 Log

# canonical ensemble production
Metropolis trials_per_cycle {tpc} cycles_to_complete {production_cycles}
Log trials_per_write {tpc} output_file {prefix}{sim}.txt
Movie trials_per_write {tpc} output_file {prefix}{sim}.xyz
Energy trials_per_write {tpc} output_file {prefix}{sim}_en.txt
CPUTime trials_per_write {tpc} output_file {prefix}{sim}_cpu.txt
Run until complete
""".format(**params))

def post_process(params):
    """ Plot energy and compare with https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm """
    ens = np.zeros(shape=(params['num_sims'], 2))
    for sim in range(params['num_sims']):
        log = pd.read_csv(params['prefix']+str(sim)+'.txt')
        assert int(log['num_particles_of_type0'][0]) == params['num_particles']
        energy = pd.read_csv(params['prefix']+str(sim)+'_en.txt')
        ens[sim] = np.array([energy['average'][0],
                             energy['block_stdev'][0]])/params['num_particles']
    # data from https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
    rhos_srsw = [0.001, 0.003, 0.005, 0.007, 0.009]
    ens_srsw = [-9.9165E-03, -2.9787E-02, -4.9771E-02, -6.9805E-02, -8.9936E-02]
    en_stds_srsw = [1.89E-05, 3.21E-05, 3.80E-05, 7.66E-05, 2.44E-05]
    plt.errorbar(rhos_srsw, ens_srsw, en_stds_srsw, fmt='+', label='SRSW')
    plt.errorbar(params['rhos'], ens[:, 0], ens[:, 1], fmt='x', label='FEASST')
    plt.xlabel(r'$\rho$', fontsize=16)
    plt.ylabel(r'$U/N$', fontsize=16)
    plt.legend(fontsize=16)
    #plt.savefig(params['prefix']+'_energy.png', bbox_inches='tight', transparent='True')
    if len(rhos_srsw) == params['num_sims']: # compare with srsw exactly
        for sim in range(params['num_sims']):
            diff = ens[sim][0] - ens_srsw[sim]
            assert np.abs(diff) < 10*np.sqrt(ens[sim][1]**2 + en_stds_srsw[sim]**2)

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
