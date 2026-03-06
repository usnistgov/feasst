"""
This is an experimental and untested method that is not currently recommended
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
    parser.add_argument('--reference_sigma', type=float, default=1.,
                        help='reference potential is a hard sphere unit diameter which is also the size of the inner hard sphere in the square well.')
    parser.add_argument('--beta', type=float, default=1./0.46, help='inverse temperature')
    parser.add_argument('--sigma_hs', type=float, default=1., help='inner hard sphere diameter')
    parser.add_argument('--num_orientations_per_pi', type=int, default=4, help='number of table orientations per 180 degrees')
    parser.add_argument('--tabsize', type=int, default=7, help='tabular size for each nested')
    parser.add_argument('--tpc', type=int, default=int(1e4), help='trials per cycle')
    parser.add_argument('--equilibration_cycles', type=int, default=int(1e2),
        help='number of cycles for equilibration')
    parser.add_argument('--production_cycles', type=int, default=int(1e2),
        help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=0.1, help='hours until termination')
    parser.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
    parser.add_argument('--cutoff', type=float, default=3, help='potential cutoff distance')
    parser.add_argument('--run_type', '-r', type=int, default=0,
        help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--seed', type=int, default=-1,
    #parser.add_argument('--seed', type=int, default=1234,
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
    params['prefix'] = 'tr6d'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    params['hard_limit_u'] = 4*(np.power(params['sigma_hs'], -12) \
                               -np.power(params['sigma_hs'], -6)) + 1 # WCA
    params['cut_wca'] = np.power(2, 1./6.)
    return params, args

# write fst script to run a simulation
def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# Mayer-sampling simulation generates training data for each pair-wise interaction
MonteCarlo
RandomMT19937 seed={seed}
#particle_type=t1:/feasst/plugin/aniso/tutorial/aniso_trimer5.txt,t2:/feasst/particle/hard_sphere.txt add_num_t1_particles=1 add_num_t2_particles=1
Configuration cubic_side_length=10 periodic=false,false,false \
    particle_type=t1:/feasst/plugin/aniso/tutorial/aniso_trimer5.txt,t2:/feasst/particle/hard_sphere.txt add_num_t1_particles=1 add_num_t2_particles=1 \
    group=fixed,mobile fixed_particle_index=0 mobile_particle_index=1
#    particle_type=fluid:/feasst/plugin/aniso/tutorial/aniso_trimer.txt add_num_fluid_particles=2 \
#    cutoff={sigma_hs} sigma={sigma_hs}
#    cutoff={cut_wca}
#cutoff=3 cutoffA_R=1.122462048309373 cutoffR_R=1.122462048309373 # anisotropic trimer
WriteModelParams output_file={prefix}{sim:03d}_model_params.csv
Record save_positions={prefix}{sim:03d}_init.xyz
Potential Model=HardSphere energy_cutoff=1e50
#Potential Model=LennardJonesCutShift
#HWH RefPotential ref=hs Model=HardSphere sigma=0 sigmaA={reference_sigma} cutoff=0 cutoffA={reference_sigma}
#RefPotential ref=hs Model=HardSphere sigma=0 sigmaA={reference_sigma} cutoff=0 cutoffA={reference_sigma} sigmaA1={reference_sigma} cutoffA1={reference_sigma}
#RefPotential ref=hs Model=HardSphere sigmaR=0 cutoffR=0 epsilonR=0
#RefPotential ref=hs Model=HardSphere sigmaR1=0 cutoffR1=0
RefPotential ref=hs Model=HardSphere sigma=0 sigmaA_H={reference_sigma}
#RefPotential ref=hs Model=HardSphere sigmaR1=0 cutoffR1=0 sigmaR=0 cutoffR=0
ThermoParams beta={beta}

# initialize MayerSampling
MayerSampling trials_per_cycle={tpc} cycles_to_complete={equilibration_cycles}
TrialTranslate new_only=true ref=hs tunable_param=1 group=mobile
#TrialRotate new_only=true ref=hs tunable_param=40 group=mobile
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# tune trial parameters
Let [write]=trials_per_write={tpc} output_file={prefix}{sim:03d}
Log [write]_mayer_eq.csv
Movie [write]_mayer_eq.xyz
CriteriaWriter [write]_b2_eq.csv
Tune
Run until=complete
Remove name=Tune,CriteriaWriter,Log,Movie

# production
CriteriaWriter [write]_b2.csv
Log [write]_mayer.csv
Movie [write]_mayer.xyz
MayerSampling trials_per_cycle={tpc} cycles_to_complete={production_cycles} training_file={prefix}{sim:03d}_training.csv
Run until=complete

# Build and output ModelRecursiveTable
# start with hard particle interactions to define the 5D rh function "contact table"
# define rh based on a high energy (hard_limit_u), HS is problematic for sampling
# First, use Rotator to generate rh according to original size
# Then, MayerSampling data is used to find the element which contains the distance furthest from predictions
# A nested table is added at that element, and then the process is repeated until criteria is met
# compare the B22 of the hard particle
# Then, use MayerSampling data of the full potential to build the energy table based on rh
Record load_positions={prefix}{sim:03d}_init.xyz
#BuildRecursiveTable mayer_training_file={prefix}{sim:03d}_training.csv hard_limit_u={hard_limit_u} size={tabsize} num_orientations_per_pi={num_orientations_per_pi} beta=1 min_criteria=0.22 output_file={prefix}{sim:03d}_table.txt verbose_file={prefix}{sim:03d}_fitting.csv
BuildRecursiveTable mayer_training_file={prefix}{sim:03d}_training.csv hard_limit_u={hard_limit_u} size={tabsize} num_orientations_per_pi={num_orientations_per_pi} beta=1 min_criteria=0.1 output_file={prefix}{sim:03d}_table.txt
#cutoff=4.154700538
# cutoff is manually set to 3 + 2L

# Test by comparing MayerSampling results to the HardSphere Potential
MonteCarlo
RandomMT19937 seed={seed}
Configuration particle_type=t1:/feasst/plugin/aniso/particle/aniso_tabular.txt,t2:/feasst/particle/hard_sphere.txt \
    cubic_side_length=10 add_num_t1_particles=1 add_num_t2_particles=1 group=second second_particle_index=1
Potential Model=TwoBodyTable VisitModelInner=VisitModelInnerRecursiveTable input_file=0_H:{prefix}{sim:03d}_table.txt ignore_energy=true
RefPotential ref=hs Model=HardSphere sigma=0 sigma0_H={reference_sigma} cutoff=0 cutoff0_H={reference_sigma}
ThermoParams beta={beta}
Metropolis

TrialTranslate new_only=true ref=hs tunable_param=1 group=second
TrialRotate new_only=true ref=hs tunable_param=40
#TrialParticlePivot new_only=true ref=hs tunable_param=40 particle_type=fluid
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# tune trial parameters
MayerSampling trials_per_cycle={tpc} cycles_to_complete={equilibration_cycles}
#Let [write]=trials_per_write={tpc} output_file={prefix}{sim:03d}
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
    ax.plot(df[0], 5*df[5], color=color, label='criteria', linestyle='dashed')
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

    #en = pd.read_csv(params['prefix']+'000_two.csv')
    #print('en', en['energy'])
    #assert np.abs(en['energy'][0] + 0.93953) < 1e-5
    b2 = read_b2(params['prefix']+'000_b2.csv')
    print('b2:', b2)
    assert np.abs(b2['second_virial_ratio'] - 2.2222) < 0.02
    b2tab = read_b2(params['prefix']+'000_b22.csv')
    print('b2 table:', b2tab)
    assert np.abs(b2tab['second_virial_ratio'] - 2.2222) < 0.05

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
