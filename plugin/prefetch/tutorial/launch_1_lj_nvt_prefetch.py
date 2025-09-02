"""
Prefetching Canonical Ensemble Monte Carlo Simulation of Lennard Jones
"""

import os
import argparse
from pyfeasst import fstio

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/', help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/lj_new.txt', help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=1./0.88, help='inverse temperature')
    parser.add_argument('--num_particles', type=int, default=500, help='number of particles')
    parser.add_argument('--density', type=float, default=0.85, help='highest number density')
    parser.add_argument('--tpc', type=int, default=int(1e4), help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=int(1e2), help='number of cycles for equilibraiton')
    parser.add_argument('--production', type=int, default=int(1e3), help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
    parser.add_argument('--procs_per_sim', type=int, default=4, help='number of processors')
    parser.add_argument('--procs_per_node', type=int, default=32, help='number of processors')
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
    params['prefix'] = 'nvt'
    params['script'] = __file__
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.98*params['hours_terminate'] - 0.0333 # terminate before queue
    assert params['procs_per_node'] % params['procs_per_sim'] == 0
    params['num_sims'] = int(params['num_nodes']*params['procs_per_node']/params['procs_per_sim'])
    params['acceptances'] = 4*[0.25, 0.5]
    assert len(params['acceptances']) == params['num_sims']
    params['cubic_side_length'] = (params['num_particles']/params['density'])**(1./3.)
    return params, args

def sim_node_dependent_params(params):
    """ Set parameters that depend upon the sim or node here. """
    params['acceptance'] = fstio.prefetch_acceptance(
        acceptance=params['acceptances'][params['sim']],
        num_processors=params['procs_per_sim'])

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
Prefetch synchronize=true
RandomMT19937 seed={seed}
Configuration cubic_side_length={cubic_side_length} particle_type=lj:{fstprt}
Potential Model=LennardJones
Potential VisitModel=LongRangeCorrections
ThermoParams beta={beta} chemical_potential=-1
Metropolis
TrialTranslate tunable_param=0.2 tunable_target_acceptance={acceptance}
#Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# grand canonical ensemble initalization
TrialAdd particle_type=lj
Run until_num_particles={num_particles}
Remove name=TrialAdd

# canonical ensemble equilibration
Metropolis trials_per_cycle={tpc} cycles_to_complete={equilibration}
Tune
CheckEnergy trials_per_update={tpc} decimal_places=8
Let [write]=trials_per_write={tpc} output_file={prefix}{sim:03d}
Log [write]_eq.csv
Run until=complete
Remove name=Tune,Log

# canonical ensemble production
Metropolis trials_per_cycle={tpc} cycles_to_complete={production}
Log [write].csv
#Movie [write].xyz
For [an]:[ext]=Energy:_en.csv,CPUTime:_cpu.csv
    [an] [write][ext] append=true
EndFor
#ProfileCPU [write]_profile.csv
#GhostTrialVolume [write]_pressure.csv trials_per_update={tpc}
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
    for sim,acc in enumerate(params['acceptances']):
        time = pd.read_csv("{}{:03d}_cpu.csv".format(params['prefix'], sim), sep='\s+', header=None)
        en = pd.read_csv("{}{:03d}_en.csv".format(params['prefix'], sim), header=None, comment="a")
        equil=50
        logt = np.log(time[1][equil:])
        logs = np.log(en[2][equil:])
        popt, pcov = curve_fit(linear_fit, logt, logs)
        df.loc[sim, 'acc'] = acc
        df.loc[sim, 'b'] = popt[0]
        color='red'
        if acc == 0.25:
            color='blue'
        plt.scatter(np.log(time[1]), np.log(en[2]), color=color)
        plt.plot(logt, linear_fit(logt, popt[0]), color=color)
    grp = df.groupby('acc')
    bmin = grp.mean()['b'][0.25].min()
    df['z'] = np.exp(2*(bmin-df['b']))
    grp = df.groupby('acc')
    df = pd.DataFrame(data={'b': grp.mean()['b'], 'b_std': grp.std()['b']/2,
                            'z': grp.mean()['z'], 'z_std': grp.std()['z']/2})
    df.to_csv(params['prefix']+'_summary.csv')
    print(df['z'][0.5])
    print('Efficiency of 0.5 acceptance w.r.t. 0.25:', df['z'][0.5], 'stdev', df['z_std'][0.5])
    #plt.show()
    assert df['z'][0.5] < 1

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
