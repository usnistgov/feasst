"""
Prefetching Gibbs Ensemble Monte Carlo Simulation of Lennard Jones
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/', help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/lj.txt', help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=1./0.85, help='inverse temperature')
    parser.add_argument('--tpc', type=int, default=int(1e4), help='trials per cycle')
    parser.add_argument('--vapor_weight', type=float, default=1, help='weight of vapor nvt trial vs liquid')
    parser.add_argument('--equilibration_cycles', type=int, default=int(1e2), help='number of cycles for equilibration')
    parser.add_argument('--production_cycles', type=int, default=int(1e3), help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
    parser.add_argument('--procs_per_node', type=int, default=32, help='number of processors per node')
    parser.add_argument('--procs_per_sim', type=int, default=4, help='number of processors per simuilation')
    parser.add_argument('--run_type', '-r', type=int, default=0, help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--seed', type=int, default=-1, help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--max_restarts', type=int, default=0, help='Number of restarts in queue')
    parser.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
    parser.add_argument('--scratch', type=str, default=None, help='Optionally write scheduled job to scratch/logname/jobid.')
    parser.add_argument('--queue_flags', type=str, default="", help='extra flags for queue (e.g., for slurm, "-p queue")')
    parser.add_argument('--node', type=int, default=0, help='node ID')
    parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
    parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

    # Convert arguments into a parameter dictionary, and add argument-dependent parameters.
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    params = vars(args)
    params['script'] = __file__
    params['prefix'] = 'gibbs'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    assert params['procs_per_node'] % params['procs_per_sim'] == 0
    params['num_sims'] = int(params['num_nodes']*params['procs_per_node']/params['procs_per_sim'])
    params['equil'] = params['equilibration_cycles']*params['tpc']
    params['double_equil'] = 2*params['equil']
    params['vapor_weights'] = 4*[0.1, 10]
    return params, args

def sim_node_dependent_params(params):
    """ Set parameters that depent upon the sim or node here. """
    exp_vapor_dens = 0.012370852011749
    exp_liquid_dens = 0.7629707597552141
    num_total = 600
    params['vapor_weight'] = params['vapor_weights'][params['sim']]
    params['vapor_fraction'] = 0.2
    params['init_vapor_num'] = int(params['vapor_fraction']*num_total)
    params['init_liquid_num'] = num_total - params['init_vapor_num']
    params['init_vapor_side_length'] = (params['init_vapor_num']/exp_vapor_dens)**(1./3.)
    params['init_liquid_side_length'] = (params['init_liquid_num']/exp_liquid_dens)**(1./3.)
    if params['procs_per_sim'] == 1:
        acc = 0.25
        vacc = 0.5
    elif params['procs_per_sim'] == 2:
        acc = 0.225
        vacc = 0.45
    elif params['procs_per_sim'] == 4:
        acc = 0.2
        vacc = 0.4
    else:
        assert False
    params['acceptance'] = fstio.prefetch_acceptance(acceptance=acc,
        num_processors=params['procs_per_sim'])
    params['vol_acceptance'] = fstio.prefetch_acceptance(acceptance=vacc,
        num_processors=params['procs_per_sim'])

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
Prefetch synchronize=true load_balance=8
RandomMT19937 seed={seed}
For [config]:[len]=vapor:{init_vapor_side_length},liquid:{init_liquid_side_length}
    Configuration name=[config] cubic_side_length=[len] particle_type=fluid:{fstprt}
    Potential Model=LennardJones config=[config]
    Potential VisitModel=LongRangeCorrections config=[config]
    RefPotential VisitModel=DontVisitModel config=[config] ref=noixn
EndFor
ThermoParams beta={beta} chemical_potential=5
Metropolis
For [config]:[param]:[weight]=vapor:7:{vapor_weight},liquid:0.2:1
    TrialTranslate weight=[weight] tunable_param=[param] tunable_target_acceptance={acceptance} config=[config]
EndFor
#Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# fill both boxes with particles
Let [write]=trials_per_write={tpc} output_file={prefix}{sim:03d}
Log [write]_fill.csv
Tune
For [config]:[num]=vapor:{init_vapor_num},liquid:{init_liquid_num}
    Movie [write]_[config]_fill.xyz config=[config]
    TrialAdd particle_type=fluid config=[config]
    Run until_num_particles=[num] config=[config]
    Remove name=TrialAdd,Movie
EndFor

## gibbs initialize vapor fraction
Metropolis trials_per_cycle=1e9 cycles_to_complete=1e9
GibbsInitialize updates_density_equil={equil} updates_per_adjust={double_equil} fraction_particles_low_density {vapor_fraction}
TrialGibbsParticleTransfer weight=2 particle_type=fluid ref=noixn print_num_accepted=true configs=vapor,liquid
TrialGibbsVolumeTransfer weight=0.002 tunable_param=20. tunable_target_acceptance={vol_acceptance} ref=noixn print_num_accepted=true configs=vapor,liquid
Log [write]_eq.csv
Tune trials_per_tune=20
Run until=complete
Remove name=GibbsInitialize,Log,Tune

# gibbs ensemble production
Metropolis trials_per_cycle={tpc} cycles_to_complete={production_cycles}
Log [write].csv
For [config]=vapor,liquid
    For [analyze]:[file]=SpecificEnergy:_en.csv,SpecificVolume:_vol.csv,Density:_dens.csv
        [analyze] [write]_[config][file] config=[config] append=true
    EndFor
EndFor
#GhostTrialVolume [write]_pressure.csv trials_per_update=1e3
#ProfileCPU [write]_profile.csv
CPUTime    [write]_cpu.csv append=true
Run until=complete
""".format(**params))

def linear_fit(x, b):
    return -0.5*x + b

def post_process(params):
    """ Compute efficiency from standard deviation of energy with time """
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.optimize import curve_fit
    zrs = params['num_sims']*[0.]
    df = pd.DataFrame(data={'wt': zrs, 'b': zrs})
    for sim,wt in enumerate(params['vapor_weights']):
        time = pd.read_csv("{}{:03d}_cpu.csv".format(params['prefix'], sim), sep='\s+', header=None)
        en = pd.read_csv("{}{:03d}_liquid_en.csv".format(params['prefix'], sim), header=None, comment="a")
        equil=500
        logt = np.log(time[1][equil:])
        logs = np.log(en[2][equil:])
        popt, pcov = curve_fit(linear_fit, logt, logs)
        df.loc[sim, 'wt'] = wt
        df.loc[sim, 'b'] = popt[0]
        color='red'
        if wt == 0.1:
            color='blue'
        plt.scatter(np.log(time[1]), np.log(en[2]), color=color)
        plt.plot(logt, linear_fit(logt, popt[0]), color=color)
    grp = df.groupby('wt')
    bmin = grp.mean()['b'][0.1].min()
    df['z'] = np.exp(2*(bmin-df['b']))
    grp = df.groupby('wt')
    df = pd.DataFrame(data={'b': grp.mean()['b'], 'b_std': grp.std()['b']/2,
                            'z': grp.mean()['z'], 'z_std': grp.std()['z']/2})
    df.to_csv(params['prefix']+'_summary.csv')
    print('Need to run longer, but apparent efficiency of vapor weight 10 vs 0.1:', df['z'][10], 'stdev', df['z_std'][10])
    #plt.show()
    assert df['z'][10] < 1

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
