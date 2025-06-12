"""
Flat-histogram simulation of single-site Lennard Jones particles in the grand canonical ensemble.
Simulate the adsorption of an LJ fluid instead of a rigid porous network of WCA particles,
as described in https://doi.org/10.1021/acs.jpcb.3c00613 .
"""

import argparse
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import macrostate_distribution

def parse():
    # Parse arguments from command line or change their default values.
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/lj_new.txt',
                        help='FEASST particle definition')
    parser.add_argument('--pore', type=str, default='/feasst/plugin/confinement/particle/porel4.txt',
                        help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=1./0.8, help='inverse temperature')
    parser.add_argument('--mu', type=float, default=-1, help='chemical potential')
    parser.add_argument('--mu_init', type=float, default=10, help='initial chemical potential')
    parser.add_argument('--max_particles', type=int, default=64, help='maximum number of particles')
    parser.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
    parser.add_argument('--min_sweeps', type=int, default=2,
                        help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
    parser.add_argument('--cubic_side_length', type=float, default=9,
                        help='cubic periodic boundary length')
    parser.add_argument('--tpc', type=int, default=int(1e5), help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=1e0, help='number of cycles for equilibration')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
    parser.add_argument('--procs_per_node', type=int, default=32, help='number of processors')
    parser.add_argument('--run_type', '-r', type=int, default=0,
                        help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--seed', type=int, default=-1,
                        help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
    parser.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
    parser.add_argument('--scratch', type=str, default=None,
                        help='Optionally write scheduled job to scratch/logname/jobid.')
    parser.add_argument('--node', type=int, default=0, help='node ID')
    parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
    parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

    # Convert arguments into a parameter dictionary, and add argument-dependent parameters.
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    params = vars(args)
    params['script'] = __file__
    params['prefix'] = 'pore'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.95*params['hours_terminate'] - 0.05 # terminate FEASST before SLURM
    params['procs_per_sim'] = 1
    params['num_sims'] = params['procs_per_node']*params['num_nodes']
    params['wca'] = 2**(1./6.)
    params['dccb_cut'] = params['cubic_side_length']/int(params['cubic_side_length']/params['wca']) # maximize inside box
    params['min1'] = params['min_particles'] + 2
    params['minimums'] = [params['min_particles'], params['min1']]
    if params['num_nodes'] == 1:
        params['windows'] = macrostate_distribution.window_exponential(
            alpha=1.0, minimums=params['minimums'], maximum=params['max_particles'],
            number=params['num_sims'], overlap=1, min_size=2)
    return params, args

def sim_node_dependent_params(params):
    """ Define parameters that are dependent on the sim or node. """
    params['min_particles'] = params['windows'][params['sim']][0]
    params['max_particles'] = params['windows'][params['sim']][1]
    params['sim_start'] = 0
    params['sim_end'] = params['num_sims'] - 1

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed={seed}
Configuration cubic_side_length={cubic_side_length} particle_type=fluid:{fstprt},pore:{pore} \
    cutoffLJ=2.5 cutoffP={wca} cutoffLJ_P={wca} add_num_pore_particles=1 \
    group=fluid,pore fluid_particle_type=fluid pore_particle_type=pore
NeighborCriteria maximum_distance=1.375 minimum_distance=0.9 site_type0=LJ site_type1=LJ
Potential EnergyMap=EnergyMapNeighborCriteria neighbor_index=0 Model=LennardJonesForceShift
RefPotential Model=LennardJonesForceShift VisitModel=VisitModelCell min_length={dccb_cut} ref=dccb
ThermoParams beta={beta} chemical_potential={mu_init}
Metropolis
TrialTranslate weight=1 particle_type=fluid tunable_param=0.2
Let [AVB_args]=particle_type=fluid site=LJ1 target_particle_type=fluid target_site=LJ1
TrialAVB2 weight=0.1 [AVB_args]
TrialAVB4 weight=0.1 [AVB_args]
CheckEnergy trials_per_update={tpc} decimal_places=4
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# write the pore xyz for visualization
Movie output_file={prefix}n{node}s{sim:03d}_pore.xyz group=pore clear_file=true
WriteStepper analyze_name=Movie
Remove name=Movie

# gcmc initialization and nvt equilibration
TrialAdd particle_type=fluid
Let [write]=trials_per_write={tpc} output_file={prefix}n{node}s{sim:03d}
Log [write]_eq.csv
Tune
Run until_num_particles={min_particles} particle_type=fluid
Remove name=TrialAdd
ThermoParams beta={beta} chemical_potential={mu}
Metropolis trials_per_cycle={tpc} cycles_to_complete={equilibration}
Run until=complete
Remove name=Tune,Log

# gcmc tm production
FlatHistogram Macrostate=MacrostateNumParticles particle_type=fluid width=1 max={max_particles} min={min_particles} \
  Bias=WLTM min_sweeps={min_sweeps} min_flatness=25 collect_flatness=20 min_collect_sweeps=1
TrialTransfer    weight=2   ref=dccb num_steps=4 particle_type=fluid
TrialTransferAVB weight=0.2 ref=dccb num_steps=4 [AVB_args]
Log [write].csv
Tune [write]_tune.csv multistate=true stop_after_cycle=1
#To print trajectories for each macrostate in separate files, add the following arguments to the "Movie" lines below: multistate true multistate_aggregate false
Movie [write]_eq.xyz stop_after_cycle=1 group=fluid
Movie [write].xyz start_after_cycle=1 group=fluid
Energy [write]_en.csv multistate=true start_after_cycle=1
ProfileCPU [write]_profile.csv
CriteriaWriter [write]_crit.csv
CriteriaUpdater trials_per_update=1e5
Run until=complete

# continue until all simulations on the node are complete
WriteFileAndCheck sim={sim} sim_start={sim_start} sim_end={sim_end} file_prefix={prefix}n{node}s file_suffix=_finished.txt output_file={prefix}n{node}_terminate.txt
Run until_file_exists={prefix}n{node}_terminate.txt trials_per_file_check={tpc}
""".format(**params))

def post_process(params):
    import pandas as pd
    df=macrostate_distribution.splice_files(prefix=params['prefix']+'n0s', suffix='_crit.csv', shift=False).dataframe()
    #print(df['delta_ln_prob'])
    delta_ln_prob_1_2_3 = [5.33710, 4.66512, 4.27458]
    #print(df['delta_ln_prob_stdev'])
    delta_ln_prob_std_1_2_3 = [0.002927, 0.006865, 0.005270]
    df = df[1:4]
    print(df)
    df['delta_ln_prob_prev'] = [5.33710, 4.66512, 4.27458]
    df['delta_ln_prob_std_prev'] = [0.002927, 0.006865, 0.005270]
    print(df)
    z_factor = 6
    diff = df[df.delta_ln_prob - df.delta_ln_prob_prev > z_factor*df.delta_ln_prob_stdev]
    print(diff)
    assert len(diff) == 0

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
