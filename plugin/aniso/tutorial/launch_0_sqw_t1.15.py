"""
This tutorial is similar to flat histogram tutorial 4, but for a square well
with anisotropic "texture" from random noise.
This is more of a mechanical test of the anisotropic potentials reproducing
an isotropic potential:
https://www.nist.gov/mml/csd/chemical-informatics-group/sat-tmmc-liquid-vapor-coexistence-properties-square-well-fluid
"""

import argparse
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import macrostate_distribution

def generate_table(num_orientations_per_pi, num_z, table_file, gamma=1., delta=0.5):
    """ generate the tabular potential """
    dt = np.pi/num_orientations_per_pi
    assert num_z > 1
    dz = 1./(num_z - 1)
    num_orientations = 0
    num_elements = 0
    fl = open(table_file, 'w')
    fl.write('site_types 1 0\n' +
             'num_orientations_per_pi ' + str(num_orientations_per_pi) + '\n' +
             'gamma ' + str(gamma) + '\n' +
             'delta ' + str(delta) + '\n' +
             'num_z ' + str(num_z) + '\n' +
             'smoothing_distance -1\n')
    #for s1 in np.arange(0, 2*np.pi + dt/2, dt): #theta
    for s1 in np.arange(0, np.pi + dt/2, dt): # if i==j, avoid x < 0 for i/j swap symmetry
        for s2 in np.arange(0, np.pi + dt/2, dt): #phi
            for e1 in np.arange(-np.pi, np.pi + dt/2, dt):
                for e2 in np.arange(0, np.pi + dt/2, dt):
                    for e3 in np.arange(-np.pi, np.pi + dt/2, dt):
                        num_orientations += 1
                        # hard code square well for testing
                        rh = 1.
                        fl.write(str(rh) + " ")
                        for z in np.arange(0, 1. + dz/2, dz):
                            num_elements += 1
                            # hard code square well for testing
                            energy = -1. + 0.0001*random.uniform(-1., 1.)
                            fl.write(str(energy) + " ")
                        fl.write("\n")
    fl.close()
    print('num_orientations', num_orientations)
    print('num_elements', num_elements)

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/plugin/aniso/particle/aniso_tabular.fstprt',
                        help='FEASST particle definition')
    parser.add_argument('--cutoff', type=float, default=1.5, help='square well cutoff')
    parser.add_argument('--dccb_cut', type=float, default=1, help='DCCB cutoff')
    parser.add_argument('--num_orientations_per_pi', type=int, default=2,
                        help='number of orientations per 90 degrees in each table angle')
    parser.add_argument('--num_z', type=int, default=2, help='number of distances in table')
    parser.add_argument('--table_file', type=str, default='dat.txt', help='table file name')
    parser.add_argument('--beta', type=float, default=1./1.15, help='inverse temperature')
    parser.add_argument('--mu', type=float, default=-3, help='chemical potential')
    parser.add_argument('--mu_init', type=float, default=10, help='initial chemical potential')
    parser.add_argument('--max_particles', type=int, default=475, help='maximum number of particles')
    parser.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
    parser.add_argument('--min_sweeps', type=int, default=2,
                        help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
    parser.add_argument('--cubic_side_length', type=float, default=9,
                        help='cubic periodic boundary length')
    parser.add_argument('--tpc', type=int, default=int(1e6), help='trials per cycle')
    parser.add_argument('--equilibration_cycles', type=int, default=0,
                        help='number of cycles for equilibration')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=0.5, help='hours until termination')
    parser.add_argument('--procs_per_node', type=int, default=32, help='number of processors')
    parser.add_argument('--run_type', '-r', type=int, default=0,
                        help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--seed', type=int, default=-1,
                        help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--max_restarts', type=int, default=1, help='Number of restarts in queue')
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
    params['prefix'] = 'sqw'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.95*params['hours_terminate'] - 0.05 # terminate FEASST before SLURM
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    generate_table(num_orientations_per_pi=args.num_orientations_per_pi,
                   num_z=args.num_z,
                   table_file=args.table_file)
    params['minimums'] = [params['min_particles']]
    if params['num_nodes'] == 1:
        params['windows'] = macrostate_distribution.window_exponential(
            alpha=2.25, minimums=params['minimums'], maximum=params['max_particles'],
            number=params['num_sims'], overlap=1, min_size=5)
    return params, args

def sim_node_dependent_params(params):
    """ Define parameters that are dependent on the sim or node. """
    params['min_particles'] = params['windows'][params['sim']][0]
    params['max_particles'] = params['windows'][params['sim']][1]
    params['sim_start'] = 0
    params['sim_end'] = params['num_sims'] - 1

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt} cutoff {cutoff}
Potential Model TwoBodyTable VisitModelInner VisitModelInnerTable table_file dat.txt
RefPotential Model HardSphere cutoff {dccb_cut} VisitModel VisitModelCell min_length {dccb_cut}
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
TrialRotate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
CheckEnergy trials_per_update {tpc} decimal_places 4
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Log trials_per_write {tpc} output_file {prefix}n{node}s{sim}_eq.csv
Tune
Run until_num_particles {min_particles}
Remove name TrialAdd
ThermoParams beta {beta} chemical_potential {mu}
Metropolis trials_per_cycle {tpc} cycles_to_complete {equilibration_cycles}
Run until complete
Remove name0 Tune name1 Log

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} \
    Bias WLTM min_sweeps {min_sweeps} min_flatness 22 collect_flatness 18 min_collect_sweeps 1
TrialTransfer weight 2 particle_type 0
Log            trials_per_write {tpc} output_file {prefix}n{node}s{sim}.csv
Movie          trials_per_write {tpc} output_file {prefix}n{node}s{sim}_eq.xyz stop_after_cycle 1
Movie          trials_per_write {tpc} output_file {prefix}n{node}s{sim}.xyz start_after_cycle 1
Tune           trials_per_write {tpc} output_file {prefix}n{node}s{sim}_tune.csv multistate true stop_after_cycle 1
Energy         trials_per_write {tpc} output_file {prefix}n{node}s{sim}_en.csv multistate true start_after_cycle 1
CriteriaWriter trials_per_write {tpc} output_file {prefix}n{node}s{sim:03d}_crit.csv
CriteriaUpdater trials_per_update 1e5
Run until complete

# continue until all simulations on the node are complete
WriteFileAndCheck sim {sim} sim_start {sim_start} sim_end {sim_end} file_prefix {prefix}n{node}s file_suffix _finished.txt output_file {prefix}n{node}_terminate.txt
Run until_file_exists {prefix}n{node}_terminate.txt trials_per_file_check {tpc}
""".format(**params))

def post_process(params):
    lnpi=macrostate_distribution.splice_files(prefix=params['prefix']+'n0s', suffix='_crit.csv', shift=False)
    lnpi.set_minimum_smoothing(25)
    #lnpi.plot(show=True)
    rw = lnpi.equilibrium()
    assert np.abs(params['beta']*params['mu'] + rw + 3.194) < 1e-2
    #lnpi.plot(show=True)
    vap, liq = lnpi.split()
    rhov = vap.average_macrostate()/params['cubic_side_length']**3
    print('rhov', rhov)
    assert np.abs(rhov - 9.723E-02) < 2e-3
    rhol = liq.average_macrostate()/params['cubic_side_length']**3
    print('rhol', rhol)
    assert np.abs(rhol - 5.384E-01) < 5e-2

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
