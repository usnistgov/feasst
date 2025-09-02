"""
Prefetching Isothermal-Isobaric Ensemble Monte Carlo Simulation of Lennard Jones
"""

import os
import argparse
import json
from pyfeasst import fstio

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/', help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/lj_new.txt', help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=1./0.88, help='inverse temperature')
    parser.add_argument('--num_particles', type=int, default=512, help='number of particles')
    parser.add_argument('--pressure', type=json.loads, default=1.237, help='pressure')
    parser.add_argument('--initial_cubic_side_length', type=int, default=8.2, help='cubic periodic boundary length')
    parser.add_argument('--tpc', type=int, default=int(1e4), help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=int(1e2), help='number of cycles for equilibraiton')
    parser.add_argument('--production', type=int, default=int(1e3), help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
    parser.add_argument('--procs_per_node', type=int, default=32, help='number of processors')
    parser.add_argument('--procs_per_sim', type=int, default=4, help='number of processors')
    parser.add_argument('--run_type', '-r', type=int, default=0, help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--seed', type=int, default=-1,help='Random number generator seed. If -1, assign random seed to each sim.')
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
    params['prefix'] = 'npt'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    assert params['procs_per_node'] % params['procs_per_sim'] == 0
    params['num_sims'] = int(params['num_nodes']*params['procs_per_node']/params['procs_per_sim'])
    params['npt_acceptances'] = 4*[0.5, 0.9]
    assert len(params['npt_acceptances']) == params['num_sims']
    return params, args

def sim_node_dependent_params(params):
    """ Set parameters that depent upon the sim or node here. """
    params['npt_acceptance'] = fstio.prefetch_acceptance(
        acceptance=params['npt_acceptances'][params['sim']],
        num_processors=params['procs_per_sim'])
    params['nvt_acceptance'] = fstio.prefetch_acceptance(
        acceptance=0.2,
        num_processors=params['procs_per_sim'])

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
Prefetch synchronize=true load_balance=8
RandomMT19937 seed={seed}
Configuration cubic_side_length={initial_cubic_side_length} particle_type=lj:{fstprt}
Potential Model=LennardJones
Potential VisitModel=LongRangeCorrections
ThermoParams beta={beta} chemical_potential=-1
Metropolis
TrialTranslate weight=1 tunable_param=0.2 tunable_target_acceptance {nvt_acceptance}
#Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}
CheckEnergy trials_per_update={tpc} decimal_places=4

# gcmc initialization
TrialAdd particle_type=lj
Let [write]=trials_per_write={tpc} output_file={prefix}{sim:03d}
Log [write]_init.csv
Tune
Run until_num_particles={num_particles}
Remove name=TrialAdd,Log,Tune

# npt equilibration
ThermoParams beta={beta} pressure={pressure}
Metropolis trials_per_cycle={tpc} cycles_to_complete={equilibration}
TrialVolume weight=0.002 tunable_param=0.2 tunable_target_acceptance={npt_acceptance}
Tune trials_per_tune=20
Log [write]_eq.csv
#Movie [write]_eq.xyz
Density [write]_density_eq.csv
Run until=complete
Remove name=Tune,Log,Density

# npt production
Metropolis trials_per_cycle={tpc} cycles_to_complete={production}
Log [write].csv
#Movie [write].xyz
For [an]:[ext]=Energy:_en.csv,CPUTime:_cpu.csv,Volume:_volume.csv
    [an] [write][ext] append=true
EndFor
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
    df = pd.DataFrame(data={'acc': zrs, 'b': zrs})
    for sim,acc in enumerate(params['npt_acceptances']):
        time = pd.read_csv("{}{:03d}_cpu.csv".format(params['prefix'], sim), sep='\s+', header=None)
        en = pd.read_csv("{}{:03d}_en.csv".format(params['prefix'], sim), header=None, comment="a")
        equil=500
        logt = np.log(time[1][equil:])
        logs = np.log(en[2][equil:])
        popt, pcov = curve_fit(linear_fit, logt, logs)
        df.loc[sim, 'acc'] = acc
        df.loc[sim, 'b'] = popt[0]
        color='red'
        if acc == 0.5:
            color='blue'
        plt.scatter(np.log(time[1]), np.log(en[2]), color=color)
        plt.plot(logt, linear_fit(logt, popt[0]), color=color)
    grp = df.groupby('acc')
    bmin = grp.mean()['b'][0.5].min()
    df['z'] = np.exp(2*(bmin-df['b']))
    grp = df.groupby('acc')
    df = pd.DataFrame(data={'b': grp.mean()['b'], 'b_std': grp.std()['b']/2,
                            'z': grp.mean()['z'], 'z_std': grp.std()['z']/2})
    df.to_csv(params['prefix']+'_summary.csv')
    print('Need to run longer, but apparent efficiency of 0.9 versus 0.5 acceptance:', df['z'][0.9], 'stdev', df['z_std'][0.9])
    #plt.show()
    assert df['z'][0.9] < 1


if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
