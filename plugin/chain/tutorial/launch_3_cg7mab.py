"""
Simulate a 7-bead mAb and compute scattering intensity.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import feasstio
from pyfeasst import scattering

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/plugin/chain/forcefield/cg7mab2.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--cubic_box_length', type=float, default=90,
                    help='cubic periodic boundary conditions')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e5),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=int(1e0),
                    help='number of iterations for equilibraiton')
PARSER.add_argument('--production_iterations', type=int, default=int(1e2),
                    help='number of iterations for production')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.2, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=1., help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=2, help='number of processors')
PARSER.add_argument('--prefix', type=str, default='cg7', help='prefix for all output file names')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--seed', type=int, default=-1,
                    help='Random number generator seed. If -1, assign random seed to each sim.')
PARSER.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
PARSER.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
PARSER.add_argument('--scratch', type=str, default=None,
                    help='Optionally write scheduled job to scratch/logname/jobid.')
PARSER.add_argument('--queue_flags', type=str, default="", help='extra flags for queue (e.g., for slurm, "-p queue")')
PARSER.add_argument('--node', type=int, default=0, help='node ID')
PARSER.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
PARSER.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

# Convert arguments into a parameter dictionary, and add argument-dependent parameters.
ARGS, UNKNOWN_ARGS = PARSER.parse_known_args()
assert len(UNKNOWN_ARGS) == 0, 'An unknown argument was included: '+str(UNKNOWN_ARGS)
PARAMS = vars(ARGS)
PARAMS['script'] = __file__
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.95*PARAMS['hours_terminate'] - 0.05 # terminate before queue
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']
PARAMS['nums'] = [30, 3]
def sim_node_dependent_params(params):
    """ Set parameters that depent upon the sim or node here. """
    params['num_particles'] = params['nums'][params['sim']]

def write_feasst_script(params, file_name):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(file_name, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt}
Potential Model HardSphere VisitModel VisitModelCell min_length 5.3283
ThermoParams beta 1 chemical_potential -1
Metropolis
TrialTranslate tunable_param 2 tunable_target_acceptance 0.2
TrialParticlePivot weight 0.5 tunable_param 0.2 tunable_target_acceptance 0.25 particle_type 0
Checkpoint file_name {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# grand canonical ensemble initalization
TrialAdd particle_type 0
Run until_num_particles {num_particles}
RemoveTrial name TrialAdd

# canonical ensemble equilibration
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Tune
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-8
Log trials_per_write {trials_per_iteration} file_name {prefix}{sim}_eq.txt
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log

# canonical ensemble production
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {production_iterations}
Log trials_per_write {trials_per_iteration} file_name {prefix}{sim}.txt
Movie trials_per_write {trials_per_iteration} file_name {prefix}{sim}.xyz
PairDistribution trials_per_update 1000 trials_per_write {trials_per_iteration} dr 0.025 file_name {prefix}{sim}_gr.csv print_intra true
Scattering trials_per_update 100 trials_per_write {trials_per_iteration} num_frequency 4 file_name {prefix}{sim}_iq.csv
Energy trials_per_write {trials_per_iteration} file_name {prefix}{sim}_en.txt
CPUTime trials_per_write {trials_per_iteration} file_name {prefix}{sim}_cpu.txt
Run until_criteria_complete true
""".format(**params))

def post_process(params):
    iq3=pd.read_csv(params['prefix']+'1_iq.csv', comment="#")
    iq30=pd.read_csv(params['prefix']+'0_iq.csv', comment="#")
    grp3 = iq3.groupby('q', as_index=False)
    grp30 = iq30.groupby('q', as_index=False)
    plt.scatter(grp3.mean()['q'], grp30.mean()['i']/grp3.mean()['i'], label='direct I(10g/L)/I(1g/L)', marker='.')
    iq3rdfft = scattering.intensity(gr_file=params['prefix']+'1_gr.csv', iq_file=params['prefix']+'1_iq.csv', num_density=3/90**3, skip=10)
    iq30rdfft = scattering.intensity(gr_file=params['prefix']+'0_gr.csv', iq_file=params['prefix']+'0_iq.csv', num_density=30/90**3, skip=10)
    plt.scatter(iq3rdfft['q'], iq30rdfft['iq']/iq3rdfft['iq'], label='rdf ft I(10g/L)/I(1g/L)')
    plt.ylabel('S', fontsize=16)
    plt.xlabel('q(1/nm)', fontsize=16)
    plt.legend()
    #plt.savefig(params['prefix']+'.png', bbox_inches='tight', transparent='True')
    assert np.abs(grp30.mean()['i'][0]/grp3.mean()['i'][0] - 0.82777) < 0.05
    assert np.abs(iq30rdfft['iq'][0]/iq3rdfft['iq'][0] - 0.7784) < 0.05

if __name__ == '__main__':
    feasstio.run_simulations(params=PARAMS,
                             sim_node_dependent_params=sim_node_dependent_params,
                             write_feasst_script=write_feasst_script,
                             post_process=post_process,
                             queue_function=feasstio.slurm_single_node,
                             args=ARGS)
