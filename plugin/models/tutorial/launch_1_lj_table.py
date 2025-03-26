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

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/lj.fstprt',
        help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=1./0.9, help='inverse temperature')
    parser.add_argument('--num_particles', type=int, default=500, help='number of particles')
    parser.add_argument('--rho_lower', type=float, default=1e-3, help='lowest number density')
    parser.add_argument('--rho_upper', type=float, default=9e-3, help='highest number density')
    parser.add_argument('--tpc', type=int, default=int(1e5), help='trials per cycle')
    parser.add_argument('--equilibration_cycles', type=int, default=int(1e1),
        help='number of cycles for equilibration')
    parser.add_argument('--production_cycles', type=int, default=int(1e3),
        help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=0.1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=0.1, help='hours until termination')
    parser.add_argument('--procs_per_node', type=int, default=5, help='number of processors')
    parser.add_argument('--table_file', type=str, default='lj_table.txt', help='table file name')
    parser.add_argument('--plot_table', type=int, default=0, help='0: no plot, 1: plot')
    parser.add_argument('--num_z', type=int, default=int(1e3), help='number of table elements')
    parser.add_argument('--inner', type=float, default=0.75, help='As described in TablePotential')
    parser.add_argument('--gamma', type=int, default=-2, help='As described in TablePotential')
    parser.add_argument('--cutoff', type=float, default=3, help='potential cutoff distance')
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
    params['prefix'] = 'lj'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    params['rhos'] = np.linspace(params['rho_lower'], params['rho_upper'], num=params['num_sims'])
    params['cubic_side_lengths'] = np.power(params['num_particles']/params['rhos'], 1./3.).tolist()
    params['rhos'] = params['rhos'].tolist()
    generate_table(params)
    if args.plot_table == 1:
        # check the energy interpolated from the table against the analytical value
        dists = np.arange(params['inner'], params['cutoff'], 0.01)
        ens = list()
        for dist in dists:
            params['displacement_test'] = dist
            run_en(params)
            df = pd.read_csv('lj.csv')
            ens.append(df['TablePotential'].values[0])
        import matplotlib.pyplot as plt
        plt.plot(dists, ens, label='table')
        plt.plot(dists, user_potential(dists), color='black', linestyle='dotted', label='analytical')
        plt.xlabel('r', fontsize=16)
        plt.ylabel('U', fontsize=16)
        plt.legend(fontsize=16)
        plt.show()
    return params, args

def sim_node_dependent_params(params):
    params['cubic_side_length'] = params['cubic_side_lengths'][params['sim']]

def user_potential(distance):
    return 4*(distance**-12 - distance**-6)

def run_en(params):
    """
    Run a feasst simulation to obtain the energy between two particles as a
    function of sepration distance (params['displacement_test'])
    """
    with open("lj_two.xyz", "w") as file1: file1.write(
"""2
-1 8 8 8
0 0 0 0
1 0 0 {displacement_test}""".format(**params))
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
""".format(**params))
    syscode = subprocess.call(params['feasst_install']+"bin/fst < launch.txt > launch.log", shell=True, executable='/bin/bash')
    if syscode > 0: sys.exit(1)

def generate_table(params):
    assert params['num_z'] > 1
    dz = 1./(params['num_z'] - 1)
    rhg = params['inner']**params['gamma']
    rcg = params['cutoff']**params['gamma']
    with open(params['table_file'], 'w') as file1:
        file1.write("""site_types 1 0\ngamma {gamma}\ninner {inner}\nnum_z {num_z}\n""".format(**params))
        for z in np.arange(0, 1 + dz/2, dz):
            if z == 0:
                distance = params['inner']
            else:
                #print(z, rcg, rhg, z*(rcg-rhg)+rhg)
                distance = (z*(rcg - rhg) + rhg)**(1./params['gamma'])
            en = user_potential(distance)
            #print('distance', distance, 'en', en)
            file1.write(str(en) + " ")

# write fst script to run a simulation
def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt}
Potential Model TablePotential table_file {table_file} VisitModel VisitModelCell
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
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
