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
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt1', type=str, default='/feasst/particle/dimer_mie_CO2.txt', help='FEASST particle definition')
    parser.add_argument('--fstprt2', type=str, default='/feasst/particle/dimer_mie_N2.txt', help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=1./258.15, help='inverse temperature (K)')
    parser.add_argument('--pressure', type=float, default=5.38, help='pressure (MPa)')
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
RandomMT19937 seed={seed}
For [config]:[length]=vapor:55,liquid:32
    Configuration name=[config] cubic_side_length=[length] particle_type=pt1:{fstprt1},pt2:{fstprt2} \
        model_param_file=/feasst/particle/mie_model_parameters.txt
    WriteModelParams output_file={prefix}{sim:03d}_[config]_model_params.txt config=[config]
    Potential Model=Mie table_size=1e4 config=[config]
    Potential VisitModel=LongRangeCorrections config=[config]
    RefPotential VisitModel=DontVisitModel config=[config] ref=noixn
EndFor
ThermoParams beta={beta} chemical_potential=5,5 pressure={pressure}
Metropolis
For [config]=vapor,liquid
    For [pt]=pt1,pt2
        For [trial]:[tunable]=Translate:0.1,ParticlePivot:0.5
            Trial[trial] weight_per_number_fraction=0.5 particle_type=[pt] tunable_param=[tunable] config=[config]
        EndFor
    EndFor
EndFor
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# fill both boxes with particles
Let [write]=trials_per_write={tpc} output_file={prefix}{sim:03d}
Log [write]_fill.csv
Tune
For [config]:[num1]:[num2]=vapor:240:260,liquid:450:50
    Movie [write]_[config]_fill.xyz config=[config]
    For [pt]:[num]=pt1:[num1],pt2:[num2]
        TrialAdd particle_type=[pt] config=[config]
        Run until_num_particles=[num] particle_type=[pt] config=[config]
        Remove name=TrialAdd
    EndFor
    Remove name=Movie
EndFor
Remove name=Tune,Log

# npt equilibrate both boxes
Metropolis trials_per_cycle={tpc} cycles_to_complete={equil_npt}
For [config]=vapor,liquid
    TrialVolume weight=0.005 tunable_param=0.2 tunable_target_acceptance=0.5 config=[config]
    Movie [write]_[config]_npt.xyz config=[config]
EndFor
Log [write]_npt.csv
ProfileCPU [write]_npt_profile.csv
CheckEnergy trials_per_update={tpc} decimal_places=6
# a new tune is required when new Trials are introduced
# decrease trials per due to infrequency of volume transfer attempts
Tune trials_per_tune=20
Run until=complete
Remove name=Tune,Log,Movie,Movie,ProfileCPU

# Gibbs equilibration
Metropolis trials_per_cycle={tpc} cycles_to_complete={equil}
For [pt]=pt1,pt2
    TrialGibbsParticleTransfer weight=0.5 particle_type=[pt] ref=noixn print_num_accepted=true configs=vapor,liquid
EndFor
Log [write]_eq.csv
For [config]=vapor,liquid
    Movie [write]_[config]_eq.xyz config=[config]
EndFor
ProfileCPU [write]_eq_profile.csv
Run until=complete
Remove name=Log,Movie,Movie,ProfileCPU

# Gibbs ensemble production
Metropolis trials_per_cycle={tpc} cycles_to_complete={production_cycles}
Log [write].csv
For [config]=vapor,liquid
    For [analyze]:[file]=Density:_dens.csv,Movie:.xyz,Energy:_en.csv,Volume:_vol.csv,ProfileCPU:_profile.csv,CPUTime:_cpu.csv
        [analyze] [write]_[config][file] config=[config]
    EndFor
    For [pt]=1,2
        NumParticles [write]_[config]_n[pt].csv particle_type=pt[pt] config=[config]
    EndFor
EndFor
Run until=complete
""".format(**params))

def compare(label, average, stdev, params, z_factor=5):
    df = pd.read_csv(params['prefix']+"000_"+label+".csv")
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
    compare("vapor_dens", rhov_rhol[0][0], rhov_rhol[1][0], params)
    compare("liquid_dens", rhov_rhol[0][1], rhov_rhol[1][1], params)
    numvapn0 = pd.read_csv(params['prefix']+'000_vapor_n1.csv')['average'][0]
    numvapn1 = pd.read_csv(params['prefix']+'000_vapor_n2.csv')['average'][0]
    numliqn0 = pd.read_csv(params['prefix']+'000_liquid_n1.csv')['average'][0]
    numliqn1 = pd.read_csv(params['prefix']+'000_liquid_n2.csv')['average'][0]
    yco2 = numvapn0/(numvapn0+numvapn1)
    print('mol frac CO2 in vapor', yco2)
    assert np.abs(yco2 - 0.536) < 0.075
    xco2 = numliqn0/(numliqn0+numliqn1)
    print('mol frac CO2 in liquid', xco2)
    assert np.abs(xco2 - 0.938) < 0.075
    eq = pd.read_csv(params['prefix']+'000_eq.csv')
    prod = pd.read_csv(params['prefix']+'000.csv')
    for conf in ['vapor', 'liquid']:
        #plt.plot(eq['volume_'+str(conf)])
        n0eq = eq['num_particles_pt1_'+str(conf)]
        n1eq = eq['num_particles_pt2_'+str(conf)]
        plt.plot(eq['trial'], n0eq/(n0eq+n1eq))
        n0 = prod['num_particles_pt1_'+str(conf)]
        n1 = prod['num_particles_pt2_'+str(conf)]
        plt.plot(prod['trial'], n0/(n0+n1))
    #plt.axhline(0.8977, linestyle='dotted', color='black')
    #plt.axhline(0.4788, linestyle='dotted', color='black')
    plt.title('MIE CO2 N2 mixture P='+str(params['pressure'])+' MPa T='+str(1./params['beta']))
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
