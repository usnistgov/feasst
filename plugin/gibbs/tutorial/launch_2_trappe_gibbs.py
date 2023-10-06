"""
Example Gibbs ensemble Monte Carlo simulation of TraPPE alkanes using FEASST.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import physical_constants
from pyfeasst import macrostate_distribution

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/particle/n-butane.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--temperature', type=float, default=350, help='temperature in Kelvin')
PARSER.add_argument('--cutoff', type=float, default=12, help='real space cutoff distance')
PARSER.add_argument('--cubic_side_length_vapor', type=float, default=63,
                    help='initial PBC of vapor')
PARSER.add_argument('--cubic_side_length_liquid', type=float, default=44.9,
                    help='initial PBC of liquid')
PARSER.add_argument('--num_particles_vapor', type=int, default=61,
                    help='initial number of particles in the vapor')
PARSER.add_argument('--num_particles_liquid', type=int, default=451,
                    help='initial number of particles in the liquid')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e5),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=int(1e1),
                    help='number of iterations for equilibration')
PARSER.add_argument('--production_iterations', type=int, default=int(1e2),
                    help='number of iterations for production')
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
PARAMS['mu_init']=10
if 'n-butane' in PARAMS['fstprt']:
    PARAMS['num_sites'] = 4
    PARAMS['molecular_weight'] = 58.12
else:
    assert False, "input new num_sites and molecular_weight into PARAMS"
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
def write_grow_file(filename, params,
                    gce, # 0: canonical moves, 1: add only for box fill, 2: gibbs transfer
                    conf, conf2=-1): # the second conf is for gibbs transfer only
    params['conf'] = conf
    params['conf2'] = conf2
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
                    bond = """bond true mobile_site {site0} anchor_site {site1} num_steps 1 reference_index 0\n""".format(**params)
                    angle = """angle true mobile_site {site0} anchor_site {site1} anchor_site2 {site2} num_steps 1 reference_index 0\n""".format(**params)
                    dihedral = """dihedral true mobile_site {site0} anchor_site {site1} anchor_site2 {site2} anchor_site3 {site3} num_steps 1 reference_index 0\n""".format(**params)

                    # full regrowth insertion/deletion
                    if trial_type == 1 and (gce == 1 or gce == 2):
                        if site == 0:
                            if gce == 2:
                                f.write("""particle_type 0 configuration_index {conf} configuration_index2 {conf2} weight 1 gibbs_transfer true site {site0} num_steps 1 reference_index 0\n""".format(**params))
                            elif gce == 1:
                                f.write("""particle_type 0 configuration_index {conf} weight 1 add true site {site0} num_steps 1 reference_index 0\n""".format(**params))
                        elif site == 1:
                            f.write(bond)
                        elif site == 2:
                            f.write(angle)
                        else:
                            f.write(dihedral)

                    # reptation
                    elif trial_type == 0 and gce == 0:
                        if site == params['num_sites'] - 1:
                            write_partial(f, bond, angle, dihedral, params)
                        else:
                            if site == 0:
                                f.write("""particle_type 0 configuration_index {conf} weight 0.25 """.format(**params))
                            f.write("""reptate true mobile_site {site0} anchor_site {site1} num_steps 1 reference_index 0\n""".format(**params))

                    # partial regrow of the last site
                    if trial_type == 2 and gce == 0:
                        if site == 0:
                            f.write("""particle_type 0 configuration_index {conf} weight 0.25 """.format(**params))
                            write_partial(f, bond, angle, dihedral, params)

                f.write("\n")

write_grow_file(filename="trappe_c0_grow_canonical.txt", params=PARAMS, gce=0, conf=0)
write_grow_file(filename="trappe_c0_grow_add.txt", params=PARAMS, gce=1, conf=0)
write_grow_file(filename="trappe_c1_grow_canonical.txt", params=PARAMS, gce=0, conf=1)
write_grow_file(filename="trappe_c1_grow_add.txt", params=PARAMS, gce=1, conf=1)
write_grow_file(filename="trappe_grow_gibbs.txt", params=PARAMS, gce=2, conf=0, conf2=1)

def write_feasst_script(params, file_name):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(file_name, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length_vapor} particle_type0 {fstprt} cutoff {cutoff}
Configuration cubic_side_length {cubic_side_length_liquid} particle_type0 {fstprt} cutoff {cutoff}
Potential Model LennardJones configuration_index 0
Potential Model LennardJones configuration_index 1
Potential Model LennardJones VisitModel VisitModelIntra intra_cut 4 configuration_index 0
Potential Model LennardJones VisitModel VisitModelIntra intra_cut 4 configuration_index 1
Potential VisitModel LongRangeCorrections configuration_index 0
Potential VisitModel LongRangeCorrections configuration_index 1
RefPotential VisitModel DontVisitModel reference_index 0 configuration_index 0
RefPotential VisitModel DontVisitModel reference_index 0 configuration_index 1
#RefPotential Model LennardJones reference_index 0 configuration_index 0
#RefPotential Model LennardJones reference_index 0 configuration_index 1
#RefPotential Model LennardJones VisitModel VisitModelIntra intra_cut 4 reference_index 0 configuration_index 0
#RefPotential Model LennardJones VisitModel VisitModelIntra intra_cut 4 reference_index 0 configuration_index 1
ThermoParams beta {beta} chemical_potential 10
Metropolis
TrialTranslate weight 0.5 tunable_param 30 configuration_index 0
TrialTranslate weight 0.5 tunable_param 1 configuration_index 1
TrialParticlePivot weight 0.25 particle_type 0 tunable_param 180 pivot_site 0 configuration_index 0
TrialParticlePivot weight 0.25 particle_type 0 tunable_param 0.4 pivot_site 0 configuration_index 1
TrialParticlePivot weight 0.25 particle_type 0 tunable_param 180 0.25 pivot_site {last_site} configuration_index 0
TrialParticlePivot weight 0.25 particle_type 0 tunable_param 0.4 pivot_site {last_site} configuration_index 1
TrialGrowFile file_name trappe_c0_grow_canonical.txt
TrialGrowFile file_name trappe_c1_grow_canonical.txt
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-4

# gcmc initialization and nvt equilibration
Movie trials_per_write {trials_per_iteration} file_name {prefix}{sim}_c0_fill.xyz configuration_index 0
Movie trials_per_write {trials_per_iteration} file_name {prefix}{sim}_c1_fill.xyz configuration_index 1
Log trials_per_write {trials_per_iteration} file_name {prefix}{sim}_fill.csv include_bonds true
Tune

# fill the first box
TrialGrowFile file_name trappe_c0_grow_add.txt
Run until_num_particles {num_particles_vapor} configuration_index 0
RemoveTrial name_contains add

# fill the second box
TrialGrowFile file_name trappe_c1_grow_add.txt
Run until_num_particles {num_particles_liquid} configuration_index 1
RemoveTrial name_contains add

# equilibrate both
#Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
#Run until_criteria_complete true
RemoveAnalyze name Log
RemoveAnalyze name Movie
RemoveAnalyze name Movie

# gibbs ensemble equilibration
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
TrialGrowFile file_name trappe_grow_gibbs.txt
TrialGibbsVolumeTransfer weight 0.001 tunable_param 3000 reference_index 0
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-8
CheckConstantVolume trials_per_update {trials_per_iteration} tolerance 1e-4
Log trials_per_write {trials_per_iteration} file_name {prefix}{sim}_eq.csv
Movie trials_per_write {trials_per_iteration} file_name {prefix}{sim}_c0_eq.xyz configuration_index 0
Movie trials_per_write {trials_per_iteration} file_name {prefix}{sim}_c1_eq.xyz configuration_index 1
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log
RemoveAnalyze name Movie
RemoveAnalyze name Movie

# gibbs ensemble production
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {production_iterations}
Log trials_per_write {trials_per_iteration} file_name {prefix}{sim}.csv
Movie trials_per_write {trials_per_iteration} file_name {prefix}{sim}_c0.xyz configuration_index 0
Movie trials_per_write {trials_per_iteration} file_name {prefix}{sim}_c1.xyz configuration_index 1
Energy trials_per_write {trials_per_iteration} file_name {prefix}{sim}_c0_en.csv configuration_index 0
Energy trials_per_write {trials_per_iteration} file_name {prefix}{sim}_c1_en.csv configuration_index 1
Density trials_per_write {trials_per_iteration} file_name {prefix}{sim}_c0_dens.csv configuration_index 0
Density trials_per_write {trials_per_iteration} file_name {prefix}{sim}_c1_dens.csv configuration_index 1
PressureFromTestVolume trials_per_update 1e3 trials_per_write {trials_per_iteration} file_name {prefix}{sim}_pressure.csv
CPUTime trials_per_write {trials_per_iteration} file_name {prefix}{sim}_cpu.csv
Run until_criteria_complete true
""".format(**params))

def post_process(params):
    z_factor = 13
    na = physical_constants.AvogadroConstant().value()
    dens_conv = 1./na*params['molecular_weight']/1e3*1e30 # convert from N/V units of molecules/A^3 to kg/m
    vapor_density = pd.read_csv(params['prefix']+"0_c0_dens.csv")
    vapor_density['average'] *= dens_conv
    vapor_density['block_stdev'] *= dens_conv
    vapor_density['diff'] = np.abs(vapor_density['average']-30.6)
    vapor_density['tol'] = np.sqrt(vapor_density['block_stdev']**2+(2**2))
    print(vapor_density)
    diverged = vapor_density[vapor_density['diff'] > z_factor*vapor_density['tol']]
    print(diverged)
    assert len(diverged) == 0
    liquid_density = pd.read_csv(params['prefix']+"0_c1_dens.csv")
    liquid_density['average'] *= dens_conv
    liquid_density['block_stdev'] *= dens_conv
    liquid_density['diff'] = np.abs(liquid_density['average']-508)
    liquid_density['tol'] = np.sqrt(liquid_density['block_stdev']**2+(2**2))
    print(liquid_density)
    diverged = liquid_density[liquid_density['diff'] > z_factor*liquid_density['tol']]
    print(diverged)
    assert len(diverged) == 0


if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
