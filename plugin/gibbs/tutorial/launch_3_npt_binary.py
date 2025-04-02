"""
Gibbs ensemble simulation of a binary mixture of CO2 and N2 using SAFT-based MIE models.
https://doi.org/10.1021/acs.jpcb.5c00536
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import physical_constants

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default=os.path.expanduser('~')+'/feasst/build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt0', type=str, default='/feasst/particle/dimer_mie_CO2.fstprt', help='FEASST particle definition')
    parser.add_argument('--fstprt1', type=str, default='/feasst/particle/dimer_mie_N2.fstprt', help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=1./258.15, help='inverse temperature (K)')
    parser.add_argument('--pressure', type=float, default=7.6148, help='pressure (MPa)')
    parser.add_argument('--tpc', type=int, default=int(1e5), help='trials per cycle')
    parser.add_argument('--equil_npt', type=int, default=int(1e1), help='number of cycles for np equilibraiton')
    parser.add_argument('--equil', type=int, default=int(1e1), help='number of cycles for Gibbs equilibraiton')
    parser.add_argument('--production_cycles', type=int, default=int(1e2),
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
    params['prefix'] = 'binary'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    # convert pressure in MPa to K/Ang^3
    # MPa(1e6 Pa/MPa)(J/m^3/Pa)(m^3/1e30/A^3)(K/kB/J)
    params['pressure'] = params['pressure']/1e24/physical_constants.BoltzmannConstant().value()
    return params, args

def sim_node_dependent_params(params):
    """ Set parameters that depent upon the sim or node here. """

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length 55 particle_type0 {fstprt0} particle_type1 {fstprt1} sigma0_1 2.9216 epsilon0_1 121.5 mie_lambda_r0_1 20.27 mie_lambda_a0_1 5.48294
Configuration cubic_side_length 32 particle_type0 {fstprt0} particle_type1 {fstprt1} sigma0_1 2.9216 epsilon0_1 121.5 mie_lambda_r0_1 20.27 mie_lambda_a0_1 5.48294
CopyFollowingLines for_num_configurations 2
    Potential Model Mie table_size 1e4
    Potential VisitModel LongRangeCorrections
    RefPotential VisitModel DontVisitModel
EndCopy
ThermoParams beta {beta} chemical_potential0 5 chemical_potential1 5 pressure {pressure}
Metropolis
CopyFollowingLines for_num_configurations 2
    TrialTranslate weight_per_number_fraction 0.5 particle_type 0 tunable_param 0.1
    TrialTranslate weight_per_number_fraction 0.5 particle_type 1 tunable_param 0.1
    TrialParticlePivot weight_per_number_fraction 0.5 particle_type 0 tunable_param 0.5
    TrialParticlePivot weight_per_number_fraction 0.5 particle_type 1 tunable_param 0.5
EndCopy
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# fill both boxes with particles
Log trials_per_write {tpc} output_file {prefix}{sim}_fill.csv
CopyNextLine replace0 configuration_index with0 0 replace1 output_file with1 {prefix}{sim}_c0_fill.xyz
Movie trials_per_write {tpc} output_file {prefix}{sim}_c1_fill.xyz configuration_index 1
Tune
TrialAdd particle_type 0 configuration_index 0
Run until_num_particles 240 particle_type 0 configuration_index 0
Remove name TrialAdd
TrialAdd particle_type 1 configuration_index 0
Run until_num_particles 260 particle_type 1 configuration_index 0
Remove name TrialAdd
TrialAdd particle_type 0 configuration_index 1
Run until_num_particles 450 particle_type 0 configuration_index 1
Remove name TrialAdd
TrialAdd particle_type 1 configuration_index 1
Run until_num_particles 50 particle_type 1 configuration_index 1
Remove name0 Tune name1 TrialAdd name2 Log name3 Movie name4 Movie

# npt equilibrate both boxes
Metropolis trials_per_cycle {tpc} cycles_to_complete {equil_npt}
CopyFollowingLines for_num_configurations 2 replace_with_index [config]
    TrialVolume weight 0.005 tunable_param 0.2 tunable_target_acceptance 0.5
    Movie trials_per_write {tpc} output_file {prefix}{sim}_c[config]_npt.xyz
EndCopy
Log trials_per_write {tpc} output_file {prefix}{sim}_npt.csv
ProfileCPU trials_per_write {tpc} output_file {prefix}{sim}_npt_profile.csv
CheckEnergy trials_per_update {tpc} decimal_places 8
# a new tune is required when new Trials are introduced
# decrease trials per due to infrequency of volume transfer attempts
Tune trials_per_tune 20
Run until complete
Remove name0 Tune name1 Log name2 Movie name3 Movie name4 ProfileCPU

# Gibbs equilibration
Metropolis trials_per_cycle {tpc} cycles_to_complete {equil}
TrialGibbsParticleTransfer weight 0.5 particle_type 0 reference_index 0 print_num_accepted true
TrialGibbsParticleTransfer weight 0.5 particle_type 1 reference_index 0 print_num_accepted true
Log trials_per_write {tpc} output_file {prefix}{sim}_eq.csv
CopyNextLine replace0 configuration_index with0 0 replace1 output_file with1 {prefix}{sim}_c0_eq.xyz
Movie trials_per_write {tpc} output_file {prefix}{sim}_c1_eq.xyz configuration_index 1
ProfileCPU trials_per_write {tpc} output_file {prefix}{sim}_eq_profile.csv
Run until complete
Remove name0 Log name1 Movie name2 Movie name3 ProfileCPU

# Bibbs ensemble production
Metropolis trials_per_cycle {tpc} cycles_to_complete {production_cycles}
Log trials_per_write {tpc} output_file {prefix}{sim}.csv
CopyFollowingLines for_num_configurations 2 replace_with_index [config]
    Density      trials_per_write {tpc} output_file {prefix}{sim}_c[config]_dens.csv
    Movie        trials_per_write {tpc} output_file {prefix}{sim}_c[config].xyz
    Energy       trials_per_write {tpc} output_file {prefix}{sim}_c[config]_en.csv
    Volume       trials_per_write {tpc} output_file {prefix}{sim}_c[config]_vol.csv
    NumParticles trials_per_write {tpc} output_file {prefix}{sim}_c[config]_n0.csv particle_type 0
    NumParticles trials_per_write {tpc} output_file {prefix}{sim}_c[config]_n1.csv particle_type 1
EndCopy
#GhostTrialVolume trials_per_write {tpc} output_file {prefix}{sim}_pressure.csv trials_per_update 1e3
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
    z_factor = 12
    rhov_rhol = [[0.00294, 0.0131], [0.0005, 0.005]] # [[avs], [stdevs]]
    compare("c0_dens", rhov_rhol[0][0], rhov_rhol[1][0], params)
    compare("c1_dens", rhov_rhol[0][1], rhov_rhol[1][1], params)
    numc0n0 = pd.read_csv(params['prefix']+'0_c0_n0.csv')['average'][0]
    numc0n1 = pd.read_csv(params['prefix']+'0_c0_n1.csv')['average'][0]
    numc1n0 = pd.read_csv(params['prefix']+'0_c1_n0.csv')['average'][0]
    numc1n1 = pd.read_csv(params['prefix']+'0_c1_n1.csv')['average'][0]
    yco2 = numc0n0/(numc0n0+numc0n1)
    print('mol frac CO2 in vapor', yco2)
    assert np.abs(yco2 - 0.46453) < 0.005
    xco2 = numc1n0/(numc1n0+numc1n1)
    print('mol frac CO2 in liquid', xco2)
    assert np.abs(xco2 - 0.89715) < 0.005
    eq = pd.read_csv(params['prefix']+'0_eq.csv')
    prod = pd.read_csv(params['prefix']+'0.csv')
    for conf in range(2):
        #plt.plot(eq['volume_config'+str(conf)])
        n0eq = eq['num_particles_of_type0_config'+str(conf)]
        n1eq = eq['num_particles_of_type1_config'+str(conf)]
        plt.plot(eq['trial'], n0eq/(n0eq+n1eq))
        n0 = prod['num_particles_of_type0_config'+str(conf)]
        n1 = prod['num_particles_of_type1_config'+str(conf)]
        plt.plot(prod['trial'], n0/(n0+n1))
    plt.axhline(0.8977, linestyle='dotted', color='black')
    plt.axhline(0.4788, linestyle='dotted', color='black')
    plt.title('MIE CO2 N2 mixture P=7.6148 MPa T=258.15')
    plt.xlabel('trials', fontsize=16)
    plt.ylabel('mole fraction CO2', fontsize=16)
    #if True:
    if False:
        plt.savefig('plot.png', bbox_inches='tight', transparent=True)

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
