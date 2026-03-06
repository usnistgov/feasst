"""
This tutorial attempts to do the same as the third tutorial for tabling trimer particles,
except that it uses MayerSampling and RecursiveTables instead.
This is an experimental and untested method that is not currently recommended
for use, because a 1 dimensional table is better off using a higher-order
interpolation strategy than linear interpolation.
Regardless, this example motivates a multi-dimensional implementation.

First, a MayerSampling simulation computes the second osmotic virial coefficient
and generates a training data set.

Then, ModelRecursiveTable is built with the training set.
To compare with https://dx.doi.org/10.1063/1.4918557,
The Boyle temperature is 0.435(5) but the table is ~0.46
"""

import subprocess
import argparse
import numpy as np
from pyfeasst import fstio

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/trimer.txt', help='FEASST particle definition')
    parser.add_argument('--reference_sigma', type=float, default=1,
                        help='reference potential is a hard sphere unit diameter which is also the size of the inner hard sphere in the square well.')
    parser.add_argument('--beta', type=float, default=1./0.47, help='inverse temperature')
    parser.add_argument('--tabsize', type=int, default=16, help='tabular size for each nested')
    parser.add_argument('--tpc', type=int, default=int(1e5), help='trials per cycle')
    parser.add_argument('--equilibration_cycles', type=int, default=int(1e1),
        help='number of cycles for equilibration')
    parser.add_argument('--production_cycles', type=int, default=int(1e1),
        help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=0.1, help='hours until termination')
    parser.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
    parser.add_argument('--cutoff', type=float, default=3, help='potential cutoff distance')
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
    params['prefix'] = 'trimer'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    return params, args

# write fst script to run a simulation
def write_feasst_script(params, script_file):
    with open(params['prefix']+"_two.xyz", "w") as file1: file1.write(
"""6
-1 10 10 10
A 0 0.577350269189626000 0
R -0.5 -0.288675134594813000 0
R  0.5 -0.288675134594813000 0
A 0 0.577350269189626000 1.1
R -0.5 -0.288675134594813000 1.1
R  0.5 -0.288675134594813000 1.1""".format(**params))
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        types = [{'type': 0, 'cutoff': 3},
                 {'type': 1, 'cutoff': np.power(2, 1./6.)}]
        for typ in types:
            params['cutoff'] = typ['cutoff']
            params['type'] = typ['type']
            myfile.write("""
# Mayer-sampling simulation generates training data for each pair-wise interaction
MonteCarlo
RandomMT19937 seed={seed}
Configuration cubic_side_length=1e10 periodic=false,false,false \
    particle_type=fluid:/feasst/particle/lj_new.txt add_num_fluid_particles=2 \
    group=second second_particle_index=1 \
    cutoff={cutoff}
Potential Model=LennardJonesCutShift
#Potential Model=HardSphere sigma=0.85
RefPotential ref=hs Model=HardSphere sigma=0 sigmaLJ={reference_sigma} cutoff=0 cutoffLJ={reference_sigma}
ThermoParams beta={beta}

# initialize MayerSampling
MayerSampling trials_per_cycle={tpc} cycles_to_complete={equilibration_cycles}
TrialTranslate new_only=true ref=hs tunable_param=1 group=second
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# tune trial parameters
Let [write{type}]=trials_per_write={tpc} output_file={prefix}{sim:03d}_{type}
Log [write{type}]_mayer_eq.csv
Movie [write{type}]_mayer_eq.xyz
CriteriaWriter [write{type}]_b2_eq.csv
Tune
Run until=complete
Remove name=Tune,CriteriaWriter,Log,Movie

# production
CriteriaWriter [write{type}]_b2.csv
Log [write{type}]_mayer.csv
Movie [write{type}]_mayer.xyz
MayerSampling trials_per_cycle={tpc} cycles_to_complete={production_cycles} training_file={prefix}{sim:03d}_{type}_training.csv
Run until=complete

# Use the training data to build a ModelRecursiveTable
BuildRecursiveTable mayer_training_file={prefix}{sim:03d}_{type}_training.csv hard_limit_u=100 size={tabsize} beta=1 min_criteria=0.05 verbose_file={prefix}{sim:03d}_{type}_fitting.csv output_file={prefix}{sim:03d}_{type}_table.txt
""".format(**params))
        myfile.write("""
# Test the ModelRecursiveTable by comparing MayerSampling results to the full potential
MonteCarlo
RandomMT19937 seed={seed}
Configuration xyz_file={prefix}_two.xyz \
    group=second second_particle_index=1 \
    particle_type=fluid:{fstprt}
#Configuration cubic_side_length=1e10 periodic=false,false,false \
#    particle_type=fluid:{fstprt} add_num_fluid_particles=2 \
#    group=second second_particle_index=1 \
#    cutoff=3 cutoffR_R=1.122462048309373 cutoffA_R=1.122462048309373
#group=second second_particle_index=1
#Potential Model=LennardJonesForceShift
#Potential Model=ModelRecursiveTable mayer_training_file=A_A:{prefix}{sim:03d}_0_training.csv,A_R:{prefix}{sim:03d}_1_training.csv,R_R:{prefix}{sim:03d}_1_training.csv hard_limit_u=100 size=16 beta=1 min_criteria=0.05 verbose_file={prefix}{sim:03d}_fitting.csv
Potential Model=ModelRecursiveTable input_file=A_A:{prefix}{sim:03d}_0_table.txt,A_R:{prefix}{sim:03d}_1_table.txt,R_R:{prefix}{sim:03d}_1_table.txt
RefPotential ref=hs Model=HardSphere sigma=0 sigmaA={reference_sigma} cutoff=0 cutoffA={reference_sigma}
ThermoParams beta={beta}
Metropolis

# obtain the energy of the initial configuration
Log output_file={prefix}{sim:03d}_two.csv max_precision=true clear_file=true
Run num_trials=1
Remove name=Log

TrialTranslate new_only=true ref=hs tunable_param=1 group=second
TrialRotate new_only=true ref=hs tunable_param=40
#TrialParticlePivot new_only=true ref=hs tunable_param=40 particle_type=fluid
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# tune trial parameters
MayerSampling trials_per_cycle={tpc} cycles_to_complete={equilibration_cycles}
Let [write]=trials_per_write={tpc} output_file={prefix}{sim:03d}
Log [write]_mayer2_eq.csv
Movie [write]_mayer2_eq.xyz
CriteriaWriter [write]_b22_eq.csv
Tune
Run until=complete
Remove name=Tune,CriteriaWriter,Log,Movie

# production
CriteriaWriter [write]_b22.csv
Log [write]_mayer2.csv
Movie [write]_mayer2.xyz
MayerSampling trials_per_cycle={tpc} cycles_to_complete={production_cycles}
Run until=complete
""".format(**params))

def pltone(ax, df, xlim, ylim, color, pts):
    ax.plot(df[0], df[3], color=color, label='table')
    ax.plot(df[0], df[5], color=color, label='criteria', linestyle='dashed')
    ax.plot(df[0], df[2], color='black', label='LJ', linestyle='dotted')
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    ax.scatter(pts, 0*pts, color=color)

def read_b2(filename):
    with open(filename) as f:
        firstline = f.readline().rstrip()
    return eval(firstline)

def post_process(params):
    import pandas as pd

    en = pd.read_csv(params['prefix']+'000_two.csv')
    print('en', en['energy'])
    assert np.abs(en['energy'][0] + 0.940225) < 1e-5
    b2 = read_b2(params['prefix']+'000_b22.csv')
    print('lj table:', b2)
    assert np.abs(b2['second_virial_ratio']) < 0.2

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 1, sharex=True)

    inner=0.75
    pts = np.arange(inner, params['cutoff'], (params['cutoff']-inner)/(params['tabsize']-1))
    #print('pts', pts)

    df1 = pd.read_csv(params['prefix']+"000_0_fitting.csv_0", header=None, delimiter=r"\s+")
    print(df1, df1[0])
    pltone(axes[0], df1, xlim=[0.9, 2], ylim=[-1, 1], color='blue', pts=pts)
    axes[0].set_ylabel(r'$U_1$', fontsize=16)

    df2 = pd.read_csv(params['prefix']+"000_0_fitting.csv_1", header=None, delimiter=r"\s+")
    pltone(axes[1], df2, xlim=[0.9, 2], ylim=[-1, 1], color='red', pts=pts)
    axes[1].set_ylabel(r'$U_2$', fontsize=16)

    df3 = pd.read_csv(params['prefix']+"000_0_fitting.csv_2", header=None, delimiter=r"\s+")
    pltone(axes[2], df3, xlim=[0.9, 2], ylim=[-1, 1], color='green', pts=pts)
    axes[2].set_ylabel(r'$U_3$', fontsize=16)

#    df4 = pd.read_csv(params['prefix']+"000_0_fitting.csv_3", header=None, delimiter=r"\s+")
#    pltone(axes[3], df4, xlim=[0.9, 2], ylim=[-1, 1], color='purple', pts=pts)
#    axes[3].set_ylabel(r'$U_4$', fontsize=16)

    plt.xlabel(r'$r$', fontsize=16)
    plt.savefig(params['prefix']+'000.png', bbox_inches='tight', transparent=True)
    #plt.show()

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
