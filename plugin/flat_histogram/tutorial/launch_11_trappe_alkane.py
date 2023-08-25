"""
Flat-histogram simulation of TraPPE alkanes in the grand canonical ensemble.
https://www.nist.gov/mml/csd/chemical-informatics-group/sat-tmmc-liquid-vapor-coexistence-properties-trappe-ua-n-butane
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import feasstio
from pyfeasst import physical_constants
from pyfeasst import macrostate_distribution

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/particle/n-butane.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--temperature', type=float, default=350, help='temperature in Kelvin')
PARSER.add_argument('--beta_mu', type=float, default=-6, help='beta time chemical potential')
PARSER.add_argument('--cutoff', type=float, default=12, help='real space cutoff distance')
PARSER.add_argument('--max_particles', type=int, default=485, help='maximum number of particles')
PARSER.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
PARSER.add_argument('--min_particles_second_window', type=int, default=30, help='minimum number of particles in the second window')
PARSER.add_argument('--min_sweeps', type=int, default=2,
                    help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
PARSER.add_argument('--cubic_side_length', type=float, default=45,
                    help='cubic periodic boundary length')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e6),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=0,
                    help='number of iterations for equilibration')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.2, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=32, help='number of processors')
PARSER.add_argument('--prefix', type=str, default='trappe', help='prefix for all output file names')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--seed', type=int, default=-1,
                    help='Random number generator seed. If -1, assign random seed to each sim.')
PARSER.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
PARSER.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
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
PARAMS['beta'] = 1./(PARAMS['temperature']*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
PARAMS['mu'] = PARAMS['beta_mu']/PARAMS['beta']
PARAMS['dccb_cut'] = 4.
PARAMS['dccb_cut'] = PARAMS['cubic_side_length']/int(PARAMS['cubic_side_length']/PARAMS['dccb_cut']) # maximize inside box
PARAMS['mu_init']=10
if 'n-butane' in PARAMS['fstprt']:
    PARAMS['num_sites'] = 4
    PARAMS['molecular_weight'] = 58.12
else:
    assert False, "input new num_sites and molecular_weight into PARMS"
PARAMS['last_site'] = PARAMS['num_sites'] - 1

def write_partial(f, bond, angle, dihedral, params):
    if params['num_sites'] == 2:
        f.write(bond)
    elif params['num_sites'] == 3:
        f.write(angle)
    elif params['num_sites'] > 3:
        f.write(dihedral)
    else:
        print('unrecognized num_sites', params['num_sites'])
        assert False

# write TrialGrowFile to include grand canonical ensemble growth and canonica ensemble reptations
def write_grow_file(filename, params, gce):
    with open(filename, 'w') as f:
        f.write("TrialGrowFile\n\n")
        for inv in [True, False]:
            for trial_type in [0, 1, 2]: # 0: reptate, 1: full regrow, 2: partial regrow
                for site in range(params['num_sites']):
                    for i in range(4):
                        sign = -1
                        if (trial_type == 0 or trial_type == 2) and site != params['num_sites'] - 1:
                            sign = 1
                        params['site'+str(i)] = site + sign*i
                        if inv:
                            params['site'+str(i)] = params['num_sites'] - site - 1 - sign*i
                    bond = """bond true mobile_site {site0} anchor_site {site1} reference_index 0 num_steps 4 reference_index 0\n""".format(**params)
                    angle = """angle true mobile_site {site0} anchor_site {site1} anchor_site2 {site2} num_steps 4 reference_index 0\n""".format(**params)
                    dihedral = """dihedral true mobile_site {site0} anchor_site {site1} anchor_site2 {site2} anchor_site3 {site3} num_steps 4 reference_index 0\n""".format(**params)

                    # full regrowth insertion/deletion
                    if trial_type == 1 and gce:
                        if site == 0:
                            f.write("""particle_type 0 weight 2 transfer true site {site0} num_steps 4 reference_index 0\n""".format(**params))
                        elif site == 1:
                            f.write(bond)
                        elif site == 2:
                            f.write(angle)
                        else:
                            f.write(dihedral)

                    # reptation
                    elif trial_type == 0 and not gce:
                        if site == params['num_sites'] - 1:
                            write_partial(f, bond, angle, dihedral, params)
                        else:
                            if site == 0:
                                f.write("""particle_type 0 weight 2 """)
                            f.write("""reptate true mobile_site {site0} anchor_site {site1} num_steps 1 reference_index 0\n""".format(**params))

                    # partial regrow of the last site
                    if not gce and trial_type == 2:
                        if site == 0:
                            f.write("""particle_type 0 weight 2 """)
                            write_partial(f, bond, angle, dihedral, params)

                f.write("\n")

write_grow_file(filename="trappe_grow_canonical.txt", params=PARAMS, gce=False)
write_grow_file(filename="trappe_grow_grand_canonical.txt", params=PARAMS, gce=True)

def write_feasst_script(params, file_name):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(file_name, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_checkpoint} ln_prob_file {prefix}n{node}_lnpi.txt min_window_size -1
WindowExponential maximum {max_particles} min0 {min_particles} min1 {min_particles_second_window} num {procs_per_node} overlap 0 alpha 2.15 min_size 3
Checkpoint file_name {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt} cutoff {cutoff}
Potential Model LennardJones
Potential Model LennardJones VisitModel VisitModelIntra intra_cut 4
Potential VisitModel LongRangeCorrections
RefPotential Model LennardJones VisitModel VisitModelCell min_length {dccb_cut} reference_index 0
RefPotential Model LennardJones VisitModel VisitModelIntra intra_cut 4 reference_index 0
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 0.25 particle_type 0 tunable_param 0.2 tunable_target_acceptance 0.25 pivot_site 0
TrialParticlePivot weight 0.25 particle_type 0 tunable_param 0.2 tunable_target_acceptance 0.25 pivot_site {last_site}
TrialGrowFile file_name trappe_grow_canonical.txt
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-4

# gcmc initialization and nvt equilibration
TrialGrowFile file_name trappe_grow_grand_canonical.txt
Log trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_eq.txt
Tune
Run until_num_particles [soft_macro_min]
RemoveTrial index 7
RemoveTrial index 7
RemoveTrial index 7
RemoveTrial index 7
ThermoParams beta {beta} chemical_potential {mu}
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} min_flatness 25 collect_flatness 20 min_collect_sweeps 1
TrialGrowFile file_name trappe_grow_grand_canonical.txt
Log trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index].txt
Movie trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_eq.xyz stop_after_iteration 1
Movie trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index].xyz start_after_iteration 1
Tune trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_tune.txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_en.txt multistate true start_after_iteration 1
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_crit.txt
""".format(**params))

def post_process(params):
    lnp = macrostate_distribution.splice_collection_matrix(prefix='trappen0s', suffix='_crit.txt', use_soft=True)
    lnp.equilibrium()
    #lnp.plot(show=True)
    print('WARNING: max_particles should be higher but the liquid peak was truncated to make the simulation faster')
    vapor, liquid = lnp.split()
    volume = params['cubic_side_length']**3
    na = physical_constants.AvogadroConstant().value()
    dens_conv = 1./volume/na*params['molecular_weight']/1e3*1e30 # convert from N/V units of molecules/A^3 to kg/m
    # https://www.nist.gov/mml/csd/chemical-informatics-research-group/sat-tmmc-liquid-vapor-coexistence-properties-trappe-ua-n
    assert np.abs(30.6 - vapor.average_macrostate()*dens_conv) < 2
    assert np.abs(508 - liquid.average_macrostate()*dens_conv) < 30

if __name__ == '__main__':
    feasstio.run_simulations(params=PARAMS,
                             sim_node_dependent_params=None,
                             write_feasst_script=write_feasst_script,
                             post_process=post_process,
                             queue_function=feasstio.slurm_single_node,
                             args=ARGS)
