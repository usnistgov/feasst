"""
Simulate hard spheres and compute scattering intensity.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import scattering

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/atom.fstprt',
                        help='FEASST particle definition')
    parser.add_argument('--num_particles', type=int, default=128, help='number of particles')
    parser.add_argument('--cubic_side_length', type=float, default=8,
                        help='cubic periodic boundary conditions')
    parser.add_argument('--tpc', type=int, default=int(1e5), help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=int(1e0), help='number of cycles for equilibration')
    parser.add_argument('--production', type=int, default=int(1e1), help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=1., help='hours until termination')
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
    params['prefix'] = 'hs'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    return params, args

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt}
Potential Model HardSphere VisitModel VisitModelCell min_length 1
ThermoParams beta 1 chemical_potential -1
Metropolis
TrialTranslate tunable_param 0.2 tunable_target_acceptance 0.25
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# grand canonical ensemble initalization
TrialAdd particle_type 0
Run until_num_particles {num_particles}
Remove name TrialAdd

# canonical ensemble equilibration
Metropolis trials_per_cycle {tpc} cycles_to_complete {equilibration}
Tune
CheckEnergy trials_per_update {tpc} tolerance 1e-8
Log trials_per_write {tpc} output_file {prefix}{sim}_eq.txt
Run until complete
Remove name0 Tune name1 Log

# canonical ensemble production
Metropolis trials_per_cycle {tpc} cycles_to_complete {production}
Log trials_per_write {tpc} output_file {prefix}{sim}.txt
Movie trials_per_write {tpc} output_file {prefix}{sim}.xyz
PairDistribution trials_per_update 1000 trials_per_write {tpc} dr 0.025 output_file {prefix}{sim}_gr.csv
Scattering trials_per_update 100 trials_per_write {tpc} num_frequency 10 output_file {prefix}{sim}_iq.csv
Run until complete
""".format(**params))

def post_process(params):
    gr=pd.read_csv(params['prefix'] + '0_gr.csv', comment="#")
    iq=pd.read_csv(params['prefix'] + '0_iq.csv', comment="#")
    grp = iq.groupby('q', as_index=False)
    assert np.abs(gr['g0-0'][45] - 1.2829) < 0.05
    assert np.abs(iq['i'][3810] - 0.0988677) < 0.4
    assert np.abs(iq['i'][0]/iq['p0'][0]**2 - 1) < 0.075

    # scale the gr closer to one at the tail by dividing by the average of the last 5%
    # A fourier transform of this scaled gr will result in a smoother low-q
    gr_scaled = gr['g0-0']/np.average(gr['g0-0'][int(0.95*len(gr['r'])):])
    plt.plot(gr['r'], gr['g0-0'], label='gr')
    plt.scatter(iq['q'], iq['i']/iq['p0']**2, label='sq_all')
    plt.plot(grp.mean()['q'], grp.mean()['i']/grp.mean()['p0']**2, label='sq_av')
    qs = np.arange(2*np.pi/8, 10, 0.01)
    sq=list(); sqs=list()
    number_density = params['num_particles']/params['cubic_side_length']**3
    for q in qs:
        sqs.append(scattering.structure_factor(gr['r'], gr_scaled, frequency=q, number_density=number_density))
        sq.append(scattering.structure_factor(gr['r'], gr['g0-0'], frequency=q, number_density=number_density))
    plt.plot(qs, sq, label='sq_from_gr')
    plt.plot(qs, sqs, label='sq_from_scaled_gr')
    #plt.savefig(params['prefix']+'.png', bbox_inches='tight', transparent='True')
    #plt.show()

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
