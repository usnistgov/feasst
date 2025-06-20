"""
Example SPC/E canonical ensemble Monte Carlo simulation using FEASST.
Approximately compare with https://doi.org/10.1063/1.476834.
There are systematic differences in the energy due to different Ewald cutoffs, etc.
"""

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
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/spce_new.txt',
        help='FEASST particle definition')
    parser.add_argument('--temperature', type=float, default=298, help='temperature in Kelvin')
    parser.add_argument('--num_particles', type=int, default=512, help='number of particles')
    parser.add_argument('--cubic_side_length', type=float, default=24.8586887,
        help='cubic periodic boundary length')
    parser.add_argument('--tpc', type=int, default=int(1e5), help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=int(1e1),
        help='number of cycles for equilibration')
    parser.add_argument('--production', type=int, default=int(1e1),
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
    params['prefix'] = 'spce'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    params['procs_per_sim'] = 1
    params['cutoff'] = 0.5*params['cubic_side_length']
    params['alpha'] = 5.6/params['cubic_side_length']
    params['beta'] = 1./(params['temperature']*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
    return params, args

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed={seed}
Configuration cubic_side_length={cubic_side_length} particle_type=spce:{fstprt} physical_constants=CODATA2010 cutoff={cutoff}
Potential VisitModel=Ewald alpha={alpha} kmax_squared=38
Potential Model=ModelTwoBodyFactory models=LennardJones,ChargeScreened VisitModel=VisitModelCutoffOuter erfc_table_size=2e4
Potential Model=ChargeScreenedIntra VisitModel=VisitModelBond
Potential Model=ChargeSelf
Potential VisitModel=LongRangeCorrections
ThermoParams beta={beta} chemical_potential=1
Metropolis
TrialTranslate weight=0.5 tunable_param=2
TrialParticlePivot weight=0.5 particle_type=spce tunable_param=0.5
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# grand canonical ensemble initalization
TrialAdd particle_type=spce
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
Movie [write].xyz
Energy [write]_en.csv
CPUTime [write]_cpu.txt
Run until=complete
""".format(**params))

def post_process(params):
    """ Approximately compare energy with https://doi.org/10.1063/1.476834 """
    if params['num_sims'] == 1: # compare energy
        log = pd.read_csv(params['prefix']+'000.csv')
        assert int(log['num_particles_spce'][0]) == params['num_particles']
        energy = pd.read_csv(params['prefix']+'000_en.csv')
        diff = energy['average'][0] - (-46.82)*params['num_particles']
        assert np.abs(diff) < 20*np.sqrt(energy['block_stdev'][0]**2 + (0.02*params['num_particles'])**2)

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
