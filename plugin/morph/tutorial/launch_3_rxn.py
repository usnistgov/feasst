"""
This tutorial is still in development.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import macrostate_distribution

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--beta', type=float, default=1, help='inverse temperature')
    parser.add_argument('--num_particles', type=int, default=20, help='total number of particles')
    parser.add_argument('--cubic_side_length', type=float, default=8,
                        help='cubic periodic boundary length')
    parser.add_argument('--tpc', type=int, default=int(1e0), help='trials per cycle')
    parser.add_argument('--equilibration_cycles', type=int, default=1e3,
                        help='number of cycles for equilibration')
    parser.add_argument('--production_cycles', type=int, default=1e3,
                        help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=0.2, help='hours until termination')
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
    params['prefix'] = 'rxn'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.95*params['hours_terminate'] - 0.05 # terminate FEASST before SLURM
    params['num_sims'] = params['num_nodes']
    params['procs_per_sim'] = params['procs_per_node']
    return params, args

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
Configuration cubic_side_length={cubic_side_length} particle_type=\
reactant1:/feasst/particle/lj_new.txt,product1:/feasst/particle/lj_new.txt,\
reactant2:/feasst/particle/lj_new.txt,product2:/feasst/particle/lj_new.txt add_num_reactant2_particles=1
Potential Model=LennardJones
Potential VisitModel=LongRangeCorrections
RefPotential VisitModel=DontVisitModel
ThermoParams beta={beta} chemical_potential=1,1,1,1
Metropolis
# Add weight per num fraction
For [pt]=reactant1,reactant2,product1,product2
    TrialTranslate particle_type=[pt] weight_per_number_fraction 0.125
#TrialParticlePivot particle_type=[pt] weight_per_number_fraction 0.125
EndFor
CheckEnergy trials_per_update={tpc} decimal_places=8
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# initialization number of particles
Let [write]=trials_per_write={tpc} output_file={prefix}n{node}s{sim:03d}
Log [write]_fill.csv
Tune
TrialAdd particle_type=reactant1
Run until_num_particles={num_particles} particle_type=reactant1
Remove name=TrialAdd,Log

# equilibration
Metropolis trials_per_cycle={tpc} cycles_to_complete={equilibration_cycles}
Let [TrialMorph]=TrialMorph weight=0.1 reference_index=0 particle_type
[TrialMorph]=reactant1,reactant2 particle_type_morph=product1,product2
[TrialMorph]=product1,product2 particle_type_morph=reactant1,reactant2
[TrialMorph]=product1,reactant1 particle_type_morph=reactant1,product1
[TrialMorph]=product2,reactant2 particle_type_morph=reactant2,product2
Log [write]_eq.csv
Run until=complete
Remove name=Tune,Log

# production
Metropolis trials_per_cycle={tpc} cycles_to_complete={production_cycles}
Log [write].csv
Tune [write]_tune.csv
Movie [write]_eq.xyz
Movie [write].xyz
Energy [write]_en.csv
ProfileCPU [write]_profile.csv trials_per_update=5e3
NumParticles [write]_n.csv particle_type=reactant2
Run until=complete
""".format(**params))

def post_process(params):
    """ Placeholder """
    assert True

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
