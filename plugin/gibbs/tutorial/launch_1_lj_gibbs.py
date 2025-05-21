"""
Example single-site Lennard-Jones Gibbs ensemble Monte Carlo simulation using FEASST.
"""

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
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/lj.txt',
                        help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=1., help='inverse temperature')
    parser.add_argument('--tpc', type=int, default=int(1e5), help='trials per cycle')
    parser.add_argument('--equilibration_cycles', type=int, default=int(1e2),
                        help='number of cycles for equilibraiton')
    parser.add_argument('--production_cycles', type=int, default=int(3e2),
                        help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
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
    params['prefix'] = 'lj'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    params['equil'] = params['equilibration_cycles']*params['tpc']
    params['double_equil'] = 2*params['equil']
    return params, args

def sim_node_dependent_params(params):
    """ Set parameters that depent upon the sim or node here. """

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
# purposefully start with a bad volume guess to see if equilibration adjusts volume
Configuration cubic_side_length 16 particle_type0 {fstprt}
Configuration cubic_side_length 8 particle_type0 {fstprt}
CopyFollowingLines for_num_configurations 2
    Potential Model LennardJones
    Potential VisitModel LongRangeCorrections
    RefPotential VisitModel DontVisitModel
EndCopy
ThermoParams beta {beta} chemical_potential 5
Metropolis
CopyNextLine replace0 configuration_index with0 0 replace1 tunable_param with1 2.0
TrialTranslate tunable_param 0.1 tunable_target_acceptance 0.2 configuration_index 1
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# fill both boxes with particles
Log trials_per_write {tpc} output_file {prefix}{sim}_fill.csv
CopyNextLine replace0 configuration_index with0 0 replace1 output_file with1 {prefix}{sim}_c0_fill.xyz
Movie trials_per_write {tpc} output_file {prefix}{sim}_c1_fill.xyz configuration_index 1
Tune
TrialAdd particle_type 0 configuration_index 0
Run until_num_particles 112 configuration_index 0
Remove name TrialAdd
TrialAdd particle_type 0 configuration_index 1
Run until_num_particles 400 configuration_index 1
Remove name0 Tune name1 TrialAdd name2 Log name3 Movie name4 Movie

# gibbs equilibration cycles: equilibrate, estimate density, adjust, repeat
# start a very long run GibbsInitialize completes once targets are reached
Metropolis trials_per_cycle 1e9 cycles_to_complete 1e9
GibbsInitialize updates_density_equil {equil} updates_per_adjust {double_equil}
TrialGibbsParticleTransfer weight 0.5 particle_type 0 reference_index 0 print_num_accepted true
TrialGibbsVolumeTransfer weight 0.01 tunable_param 10. tunable_target_acceptance 0.5 reference_index 0 print_num_accepted true
CheckEnergy trials_per_update {tpc} decimal_places 8
Log trials_per_write {tpc} output_file {prefix}{sim}_eq.csv
CopyNextLine replace0 configuration_index with0 0 replace1 output_file with1 {prefix}{sim}_c0_eq.xyz
Movie trials_per_write {tpc} output_file {prefix}{sim}_c1_eq.xyz configuration_index 1
ProfileCPU trials_per_write {tpc} output_file {prefix}{sim}_eq_profile.csv
# a new tune is required when new Trials are introduced
# decrease trials per due to infrequency of volume transfer attempts
Tune trials_per_tune 20
Run until complete
Remove name0 GibbsInitialize name1 Tune name2 Log name3 Movie name4 Movie name5 ProfileCPU

# gibbs ensemble production
Metropolis trials_per_cycle {tpc} cycles_to_complete {production_cycles}
Log trials_per_write {tpc} output_file {prefix}{sim}.csv
CopyFollowingLines for_num_configurations 2 replace_with_index [config]
    Density trials_per_write {tpc} output_file {prefix}{sim}_c[config]_dens.csv
    Movie   trials_per_write {tpc} output_file {prefix}{sim}_c[config].xyz
    Energy  trials_per_write {tpc} output_file {prefix}{sim}_c[config]_en.csv
    Volume  trials_per_write {tpc} output_file {prefix}{sim}_c[config]_vol.csv
EndCopy
GhostTrialVolume trials_per_write {tpc} output_file {prefix}{sim}_pressure.csv trials_per_update 1e3
ProfileCPU trials_per_write {tpc} output_file {prefix}{sim}_profile.csv
CPUTime    trials_per_write {tpc} output_file {prefix}{sim}_cpu.csv
Run until complete
""".format(**params))

def compare(label, average, stdev, params, z_factor=5):
    df = pd.read_csv(params['prefix']+"0_"+label+".csv")
    df['diff'] = np.abs(df['average']-average)
    df['tol'] = np.sqrt(df['block_stdev']**2+stdev**2)
    print(label, df)
    diverged = df[df['diff'] > z_factor*df['tol']]
    if len(diverged) > 0:
        print(diverged)
    assert len(diverged) == 0

def post_process(params):
    z_factor = 3
    #fh rhov_rhol_p = [[0.1003, 0.56329, 0.07721], [9.41E-06, 4.51E-05, 5.7E-06]] # T=1.2 srsw fh
    rhov_rhol_p = [[2.9556E-02, 7.0094E-01, 2.4950E-02], [3.45E-06, 6.31E-05, 1.67E-06]] # T=1 srsw fh
    #fh rhov_rhol_p = [[6.1007E-03, 0.79981, 0.0046465], [5.63E-07, 0.000013, 3.74e-7]] #T=0.8 srsw fh
    compare("c0_dens", rhov_rhol_p[0][0], rhov_rhol_p[1][0], params)
    compare("c1_dens", rhov_rhol_p[0][1], rhov_rhol_p[1][1], params)
    compare("pressure", rhov_rhol_p[0][2], rhov_rhol_p[1][2], params)
    #if True: # set to true to plot
    if False: # set to true to plot
        df = pd.read_csv('lj0_eq.csv')
        print(df)
        #plt.plot(df['volume_config0'])
        label='num_particles_of_type0'
        #label='volume'
        #label='energy'
        if label != 'num_particles_of_type0':
            for config in ['0', '1']:
                plt.plot(df[label+'_config'+config], label=config)
            plt.ylabel(label, fontsize=16)
        else:
            frac_vapor = df[label+'_config0']/(df[label+'_config0']+df[label+'_config1'])
            plt.plot(frac_vapor)
            plt.ylabel('number fraction in vapor', fontsize=16)
            plt.axhline(0.15)
            plt.axhline(0.1, linestyle='dashed')
            plt.axhline(0.2, linestyle='dashed')
        plt.xlabel('trials / 1e5', fontsize=16)
        plt.savefig('plot.png', bbox_inches='tight')

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
