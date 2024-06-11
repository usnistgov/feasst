"""
Compute the second virial coefficient using coarse-grained models.
"""

import argparse
import numpy as np
import subprocess
from pyfeasst import fstio
from pyfeasst import physical_constants

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--temperature', type=float, default=298.15, help='temperature in Kelvin')
PARSER.add_argument('--table_file', type=str, default='energy.txt', help='file describe cg table potential.')
PARSER.add_argument('--cubic_side_length', type=float, default=1e4, help='cubic side length')
PARSER.add_argument('--molecular_weight', type=float, default=14315.03534, help='molecular weight of protein in g/mol')
PARSER.add_argument('--reference_sigma', type=float, default=30, help='size of hard sphere on COM of rigid domain')
PARSER.add_argument('--ignore_energy', type=str, default='false', help='true if interaction is excluded volume only.')
PARSER.add_argument('--ignore_intra_energy', type=str, default='false', help='true if intra interaction is excluded volume only.')
PARSER.add_argument('--num_beta_taylor', type=int, default=10, help='number of Tayler series derivatives')
PARSER.add_argument('--show_plot', type=int, default=0, help='show extrapolation plot if != 0')
PARSER.add_argument('--trials_per', type=int, default=1e5, help='number of trials per iteration')
PARSER.add_argument('--equilibration', type=int, default=2e1, help='number of iterations in equilibration')
PARSER.add_argument('--production', type=int, default=2e1, help='number of iterations in production')
PARSER.add_argument('--seed', type=int, default=-1,
                    help='Random number generator seed. If -1, assign random seed to each sim.')
PARSER.add_argument('--fstprt', type=str, default='/feasst/plugin/aniso/particle/aniso_tabular.fstprt', help='fstprt file')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--prefix', type=str, default='b2', help='prefix for all output file names')
PARSER.add_argument('--hours_terminate', type=float, default=5*24, help='hours until termination')
PARSER.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
PARSER.add_argument('--procs_per_node', type=int, default=32, help='Number of nodes in queue')
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
PARAMS['sim_id_file'] = PARAMS['prefix'] + '_sim_ids.txt'
PARAMS['script'] = __file__
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.99*PARAMS['hours_terminate'] - 0.0333 # terminate before queue
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']
PARAMS['beta'] = 1./(PARAMS['temperature']*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt} add_particles_of_type0 2 group0 mobile mobile_particle_index 1
Potential Model TwoBodyTable VisitModelInner VisitModelInnerTable table_file {table_file} ignore_energy {ignore_energy}
RefPotential Model HardSphere sigma 0 cutoff 0 sigma0 {reference_sigma} cutoff0 {reference_sigma}
ThermoParams beta {beta}
MayerSampling num_beta_taylor {num_beta_taylor} num_trials_per_iteration {trials_per} num_iterations_to_complete {equilibration}
TrialTranslate new_only true reference_index 0 tunable_param 35 group mobile
TrialRotate new_only true reference_index 0 tunable_param 40

# tune trial parameters
CriteriaWriter trials_per_write {trials_per} output_file {prefix}_{sim}_b2_eq.txt
Log trials_per_write {trials_per} output_file {prefix}_{sim}_eq.txt
#Movie trials_per_write {trials_per} output_file {prefix}_{sim}_eq.xyz
#CPUTime trials_per_write {trials_per} output_file {prefix}_{sim}_cpu.txt append true
Tune
#Run num_trials {equilibration}
Run until_criteria_complete true
RemoveAnalyze name CriteriaWriter
RemoveAnalyze name Log
#RemoveAnalyze name Movie
RemoveModify name Tune

# production
CriteriaWriter trials_per_write {trials_per} output_file {prefix}_{sim}_b2.txt
Log trials_per_write {trials_per} output_file {prefix}_{sim}.txt
#Movie trials_per_write {trials_per} output_file {prefix}_{sim}.xyz
MayerSampling num_beta_taylor {num_beta_taylor} num_trials_per_iteration {trials_per} num_iterations_to_complete {production}
Run until_criteria_complete true
""".format(**params))

def post_process(params):
    import pandas as pd
    from pyfeasst import accumulator
    b2acc = accumulator.Accumulator()
    ref = (2*np.pi/3)*params['reference_sigma']**3
    mw = params['molecular_weight']/1e3 # kDa
    ref *= 1e-26*physical_constants.AvogadroConstant().value()/mw/mw # A3 to 10^4 molml/g2
    for p in range(params['procs_per_node']):
        with open(params['prefix']+'_'+str(p)+"_b2.txt", 'r') as file1:
#        with open(params['prefix']+'_'+str(p)+"_b2_eq.txt", 'r') as file1:
            lines = file1.readlines()
        if len(lines) != 0:
            #print('p', p, 'lines[0]', lines[0])
            exec('iprm=' + lines[0], globals())
            #print(iprm)
            #print('ref(A^3)', ref)
            b2_molmlg2=iprm['second_virial_ratio']*ref
            #print('b2(A^3)', b2)
            #print(b2_molmlg2, '10^-4 mol*ml*g^-2')
            b2acc.add(b2_molmlg2)
            if p == 0:
                df = pd.DataFrame(data={p: iprm['beta_taylor']})
            else:
                df[p] = iprm['beta_taylor']
    print('b2 (mol*ml/g2)', b2acc.mean(), b2acc.stdev()/np.sqrt(params['procs_per_node']))
    if params['molecular_weight'] == 14315.03534:
        assert np.abs(b2acc.mean() - 1.72) < 0.1

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          #sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
