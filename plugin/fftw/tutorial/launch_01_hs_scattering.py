"""
Test ScatteringFFTW against Scattering for a simple hard-sphere system.
"""

import argparse
from pyfeasst import fstio

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/', help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/atom.txt', help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=1, help='inverse temperature')
    parser.add_argument('--cubic_side_length', type=float, default=8, help='periodic boundary condition side length')
    parser.add_argument('--num_particles', type=int, default=128, help='number of particles')
    parser.add_argument('--tpc', type=int, default=int(1e6), help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=int(1e5), help='number of cycles for equilibraiton')
    parser.add_argument('--production', type=int, default=int(1e6), help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=5*24, help='hours until termination')
    parser.add_argument('--procs_per_node', type=int, default=2, help='number of processors')
    parser.add_argument('--run_type', '-r', type=int, default=0,
                        help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--seed', type=int, default=-1,
                        help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
    parser.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
    parser.add_argument('--scratch', type=str, default=None,
                        help='Optionally write scheduled job to scratch/logname/jobid.')
    parser.add_argument('--queue_flags', type=str, default="",
                        help='extra flags for queue (e.g., for slurm, "-p queue")')
    parser.add_argument('--node', type=int, default=0, help='node ID')
    parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
    parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

    # Convert arguments into a parameter dictionary, and add argument-dependent parameters.
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    params = vars(args)
    params['prefix'] = 'lj'
    params['script'] = __file__
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.98*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    return params, args

def sim_node_dependent_params(params):
    """ Set parameters that depend upon the sim or node here. """
    if params['sim'] == 0:
        params['num_particles'] = 1

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed={seed}
Configuration cubic_side_length={cubic_side_length} particle_type=fluid:{fstprt}
Potential Model=HardSphere VisitModel=VisitModelCell min_length=1
ThermoParams beta={beta} chemical_potential=1
Metropolis
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}
TrialTranslate tunable_param=0.2
Let [write]=trials_per_write={tpc} output_file={prefix}s{sim:03d}
Log [write]_eq.csv
Tune
CheckEnergy trials_per_update={tpc} decimal_places=6

# gcmc initialization and nvt equilibration
TrialAdd particle_type=fluid
Run until_num_particles={num_particles}
Remove name=TrialAdd
Run num_trials={equilibration}

# nvt production
Movie [write].xyz
PairDistribution [write]_gr.csv trials_per_update=1000 dr=0.025
Scattering [write]_iq.csv trials_per_update=100 num_frequency=10
ScatteringFFTW [write]_iq_fftw.csv trials_per_update=100 bin_spacing=0.1 delta_rho=1
Run num_trials={production}
""".format(**params))

def post_process(params):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from pyfeasst import scattering
    iq1=pd.read_csv(params['prefix']+'s001_iq.csv', comment="#")
    grp1 = iq1.groupby('q', as_index=False)
    assert np.abs(iq1['i'][3810] - 0.1) < 0.004
    iqfftw0 = pd.read_csv(params['prefix']+'s000_iq_fftw.csv')
    iqfftw1 = pd.read_csv(params['prefix']+'s001_iq_fftw.csv')
    assert np.abs(iqfftw1['i'][6] - 1.1) < 0.005

    plt.scatter(iq1['q'], iq1['i']/iq1['p0']**2, label='sq_all')
    plt.plot(grp1.mean()['q'], grp1.mean()['i']/grp1.mean()['p0']**2, label='sq_av', color='black')
    plt.scatter(iqfftw1['q'], iqfftw1['i']/iqfftw0['i'], label='sq fftw', color='red')
    plt.xlim([0,15])
    plt.legend()
    plt.xlabel('q', fontsize=16)
    plt.ylabel('S', fontsize=16)
    #plt.show()

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
