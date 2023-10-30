"""
This tutorial is similar to tutorial 4, but for this low temperature simulation,
we will split the simulation into two different nodes.
The first node will have less particles but a higher number of sweeps required.
The second node will have dccb but not avb.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import macrostate_distribution
from pyfeasst import multistate_accumulator

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/particle/lj.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--beta', type=float, default=1./0.7, help='inverse temperature')
PARSER.add_argument('--mu', type=float, default=-4.1603632, help='chemical potential')
PARSER.add_argument('--mu_init', type=float, default=10, help='initial chemical potential')
PARSER.add_argument('--num_particles', type=int, default=475, help='number of particles')
PARSER.add_argument('--num_particles_first_node', type=int, default=375,
                    help='number of particles in the first node')
PARSER.add_argument('--cubic_side_length', type=float, default=8,
                    help='cubic periodic boundary length')
PARSER.add_argument('--dccb_cut', type=float, default=1,
                    help='dual-cut configurational bias cutoff')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e6),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=0,
                    help='number of iterations for equilibration')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.1, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=5*24, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=32, help='number of processors')
PARSER.add_argument('--prefix', type=str, default='lj_lowt', help='prefix for output file names')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--seed', type=int, default=-1,
                    help='Random number generator seed. If -1, assign random seed to each sim.')
PARSER.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
PARSER.add_argument('--num_nodes', type=int, default=2, help='Number of nodes in queue')
PARSER.add_argument('--scratch', type=str, default=None,
                    help='Optionally write scheduled job to scratch/logname/jobid.')
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
PARAMS['hours_terminate'] = 0.95*PARAMS['hours_terminate'] - 0.05 # terminate FEASST before SLURM
PARAMS['hours_terminate'] *= PARAMS['procs_per_node'] # real time -> cpu time
PARAMS['hours_checkpoint'] *= PARAMS['procs_per_node']
PARAMS['num_sims'] = PARAMS['num_nodes']
PARAMS['procs_per_sim'] = PARAMS['procs_per_node']
PARAMS['dccb_cut'] = PARAMS['cubic_side_length']/int(PARAMS['cubic_side_length']/PARAMS['dccb_cut'])
def sim_node_dependent_params(params):
    """ Define parameters that are dependent on the sim or node. """
    if params['node'] == 0:
        params['min_particles'] = 0
        params['max_particles'] = params['num_particles_first_node']
        params['gce_trial'] = 'TrialTransfer weight 2 particle_type 0\nTrialTransferAVB weight 0.2 particle_type 0'
        params['lj_potential'] = 'Potential EnergyMap EnergyMapNeighborCriteria neighbor_index 0 Model LennardJones'
        params['ref_potential'] = ''
        params['avb_trials'] = 'TrialAVB2 weight 0.1 particle_type 0\nTrialAVB4 weight 0.1 particle_type 0'
        params['min_sweeps'] = 20
        params['window_alpha'] = 2
        params['min_window_size'] = 5
    elif params['node'] == 1:
        params['min_particles'] = params['num_particles_first_node']
        params['max_particles'] = params['num_particles']
        params['gce_trial'] = 'TrialTransfer weight 2 particle_type 0 reference_index 0 num_steps 10'
        params['lj_potential'] = 'Potential Model LennardJones'
        params['ref_potential'] = """RefPotential Model LennardJones cutoff {dccb_cut} VisitModel VisitModelCell min_length {dccb_cut}""".format(**params)
        params['avb_trials'] = ''
        params['min_sweeps'] = 2
        params['window_alpha'] = 1
        params['min_window_size'] = 3

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_checkpoint} ln_prob_file {prefix}n{node}_lnpi.txt min_window_size -1
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha {window_alpha} min_size {min_window_size}
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt}
NeighborCriteria maximum_distance 1.375 minimum_distance 0.9
{lj_potential}
{ref_potential}
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
{avb_trials}
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-4

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Log trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_eq.txt
Tune
Run until_num_particles [soft_macro_min]
RemoveTrial name TrialAdd
ThermoParams beta {beta} chemical_potential {mu}
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} min_flatness 22 collect_flatness 18 min_collect_sweeps 1
{gce_trial}
Log trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index].txt
Movie trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_eq.xyz stop_after_iteration 1
Movie trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index].xyz start_after_iteration 1
Tune trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_tune.txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_en.txt multistate true start_after_iteration 1
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_crit.txt
""".format(**params))

def post_process(params):
    """ Compare the lnpi with the SRSW. """
    num_block = 32
    beta_mu_eq = list()
    rho_vapor = list()
    rho_liquid = list()
    en_vapor = list()
    en_liquid = list()
    pressure = list()
    for block in range(-1, num_block):
        if block == -1:
            ln_prob_header = 'ln_prob'
            energy_header = 'e_average'
        else:
            ln_prob_header = 'ln_prob' + str(block)
            energy_header = 'e_block' + str(block)
        lnpi = macrostate_distribution.splice_files(prefix=params['prefix']+'n', suffix='_lnpi.txt',
                                                    ln_prob_header=ln_prob_header)
        lnpi.set_minimum_smoothing(30)
        energy = multistate_accumulator.splice_by_node(prefix=params['prefix']+'n',
                                                       suffix='_en.txt',
                                                       num_nodes=params['num_nodes'])
        lnpi.concat_dataframe(dataframe=energy, add_prefix='e_')
        delta_beta_mu = lnpi.equilibrium()
        beta_mu_eq.append(params['beta']*params['mu'] + delta_beta_mu)
        for index, phase in enumerate(lnpi.split()):
            n_gce = phase.average_macrostate()
            rho = n_gce/params['cubic_side_length']**3
            energy = phase.ensemble_average(energy_header)/n_gce
            if index == 0:
                pressure.append(-phase.ln_prob()[0]/params['beta']/params['cubic_side_length']**3)
                rho_vapor.append(rho)
                en_vapor.append(energy)
            else:
                rho_liquid.append(rho)
                en_liquid.append(energy)
    data = pd.DataFrame(data={'rho_vapor': rho_vapor, 'rho_liquid': rho_liquid,
                              'pressure': pressure, 'en_vapor': en_vapor, 'en_liquid': en_liquid,
                              'beta_mu_eq': beta_mu_eq})
    data.to_csv(params['script']+'.csv')
    for col in data.columns:
        print(col, data[col][0], '+/-', data[col][1:].std()/np.sqrt(len(data[col][1:])),
              data[col][1:].mean())

    # skip the following checks if temperature is not 0.7
    if np.abs(params['beta'] - 1./0.7) > 1e-5:
        return

    # check equilibrium properties from https://www.nist.gov/mml/csd/chemical-informatics-group/sat-tmmc-liquid-vapor-coexistence-properties-long-range
    z_factor = 6
    assert np.abs(data['rho_vapor'][0] - 1.996E-03) < z_factor*np.sqrt((data['rho_vapor'][1:].std()/np.sqrt(num_block))**2 + 1.422E-05**2)
    assert np.abs(data['rho_liquid'][0] - 8.437E-01), z_factor*np.sqrt((data['rho_liquid'][1:].std()/np.sqrt(num_block))**2 + 2.49E-04**2)
    assert np.abs(data['pressure'][0] - 1.370E-03), z_factor*np.sqrt((data['pressure'][1:].std()/np.sqrt(num_block))**2 + 9.507E-07**2)
    assert np.abs(data['en_vapor'][0] - -2.500E-02), z_factor*np.sqrt((data['en_vapor'][1:].std()/np.sqrt(num_block))**2 + 3.323E-05**2)
    assert np.abs(data['en_liquid'][0] - -6.106E+00), z_factor*np.sqrt((data['en_liquid'][1:].std()/np.sqrt(num_block))**2 + 1.832E-03**2)
    assert np.abs(data['beta_mu_eq'][0] - -6.257E+00), z_factor*np.sqrt((data['beta_mu_eq'][1:].std()/np.sqrt(num_block))**2 + 4.639E-04**2)

    # check lnpi
    lnpi = macrostate_distribution.splice_files(prefix=params['prefix']+'n', suffix='_lnpi.txt')
    srsw = pd.read_csv(params['feasst_install']+'../plugin/flat_histogram/test/data/stat070.csv')
    srsw = pd.concat([lnpi.dataframe(), srsw], axis=1)
    srsw['deltalnPI'] = srsw.lnPI-srsw.lnPI.shift(1)
    srsw.to_csv(params['prefix']+'_lnpi.csv')
    diverged = srsw[srsw.deltalnPI-srsw.delta_ln_prob > z_factor*srsw.delta_ln_prob_stdev]
    print(diverged)
    assert len(diverged) == 0

    # plot lnpi
    fst = pd.read_csv(params['prefix']+'_lnpi.csv')
    plt.plot(fst['state'], fst['ln_prob'], label='FEASST')
    plt.plot(srsw['N'], srsw['lnPI'], linestyle='dashed', label='SRSW')
    plt.xlabel('number of particles', fontsize=16)
    plt.ylabel('ln probability', fontsize=16)
    plt.legend(fontsize=16)
    #plt.savefig(params['prefix']+'_lnpi.png', bbox_inches='tight', transparent='True')

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
