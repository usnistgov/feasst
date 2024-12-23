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
PARSER.add_argument('--xyz_vapor', type=str, default='',
                    help='optionally, use an xyz file to initialize vapor')
PARSER.add_argument('--xyz_liquid', type=str, default='',
                    help='optionally, use an xyz file to initialize liquid')
PARSER.add_argument('--tpc', type=int, default=int(1e5), help='trials per cycle')
PARSER.add_argument('--equilibration_cycles', type=int, default=int(1e1),
                    help='number of cycles for equilibration')
PARSER.add_argument('--production_cycles', type=int, default=int(1e2),
                    help='number of cycles for production')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.2, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
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
PARAMS['prefix'] = 'trappe'
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.95*PARAMS['hours_terminate'] - 0.05 # terminate FEASST before SLURM
PARAMS['hours_terminate'] *= PARAMS['procs_per_node'] # real time -> cpu time
PARAMS['hours_checkpoint'] *= PARAMS['procs_per_node']
PARAMS['num_sims'] = PARAMS['num_nodes']
PARAMS['procs_per_sim'] = PARAMS['procs_per_node']
PARAMS['beta'] = 1./(PARAMS['temperature']*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
PARAMS['mu_init']=10
PARAMS['equil'] = PARAMS['equilibration_cycles']*PARAMS['tpc']
PARAMS['double_equil'] = 2*PARAMS['equil']
if 'n-butane' in PARAMS['fstprt']:
    PARAMS['num_sites'] = 4
    PARAMS['molecular_weight'] = 58.12
else:
    assert False, "input new num_sites and molecular_weight into PARAMS"
PARAMS['last_site'] = PARAMS['num_sites'] - 1
if PARAMS['xyz_vapor'] == '':
    PARAMS['vapor_config'] = """cubic_side_length {cubic_side_length_vapor}""".format(**PARAMS)
    PARAMS['init_vapor'] = """TrialGrowFile grow_file {prefix}_c0_grow_add.txt
Run until_num_particles {num_particles_vapor} configuration_index 0
Remove name_contains add""".format(**PARAMS)
else:
    PARAMS['vapor_config'] = """xyz_file {xyz_vapor}""".format(**PARAMS)
    PARAMS['init_vapor'] = ''
if PARAMS['xyz_liquid'] == '':
    PARAMS['liquid_config'] = """cubic_side_length {cubic_side_length_liquid}""".format(**PARAMS)
    PARAMS['init_liquid'] = """TrialGrowFile grow_file {prefix}_c1_grow_add.txt
Run until_num_particles {num_particles_liquid} configuration_index 1
Remove name_contains add""".format(**PARAMS)
else:
    PARAMS['liquid_config'] = """xyz_file {xyz_liquid}""".format(**PARAMS)
    PARAMS['init_liquid'] = ''

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
            for trial_type in range(3+int(params['num_sites']/2)): # 0: reptate, 1: full regrow, 2+: partial regrow
                for site in range(params['num_sites']):
                    for i in range(4):
                        sign = -1
                        if trial_type == 0 and site != params['num_sites'] - 1:
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
                                f.write("""particle_type 0 configuration_index {conf} configuration_index2 {conf2} weight 1 gibbs_transfer true site {site0} num_steps 1 reference_index 0 print_num_accepted true\n""".format(**params))
                            elif gce == 1:
                                f.write("""particle_type 0 configuration_index {conf} weight 1 add true site {site0} num_steps 1 reference_index 0\n""".format(**params))
                        elif site == 1:
                            f.write(bond)
                        elif site == 2:
                            f.write(angle)
                        else:
                            f.write(dihedral)

#                    # reptation. There seems to be a problem with reptation.
#                    elif trial_type == 0 and gce == 0:
#                        if site == params['num_sites'] - 1:
#                            write_partial(f, bond, angle, dihedral, params)
#                        else:
#                            if site == 0:
#                                f.write("""particle_type 0 configuration_index {conf} weight 0.25 """.format(**params))
#                            f.write("""reptate true mobile_site {site0} anchor_site {site1} num_steps 1 reference_index 0\n""".format(**params))
#
#                    # partial regrow of the last site
#                    if trial_type == 2 and gce == 0:
#                        if site == 0:
#                            f.write("""particle_type 0 configuration_index {conf} weight 0.25 """.format(**params))
#                            write_partial(f, bond, angle, dihedral, params)

                    # partial regrow
                    if not gce and trial_type > 1:
                        num_grow = trial_type - 1
                        if params['num_sites'] - site < num_grow:
                            if params['num_sites'] - site == num_grow - 1:
                                f.write('particle_type 0 weight '+str(2/(trial_type-2))+' ')
                            if site == 1:
                                f.write(bond)
                            elif site == 2:
                                f.write(angle)
                            elif site != 0:
                                f.write(dihedral)

                f.write("\n")

write_grow_file(filename=PARAMS['prefix']+"_c0_grow_canonical.txt", params=PARAMS, gce=0, conf=0)
if PARAMS['xyz_vapor'] == '':
    write_grow_file(filename=PARAMS['prefix']+"_c0_grow_add.txt", params=PARAMS, gce=1, conf=0)
write_grow_file(filename=PARAMS['prefix']+"_c1_grow_canonical.txt", params=PARAMS, gce=0, conf=1)
if PARAMS['xyz_liquid'] == '':
    write_grow_file(filename=PARAMS['prefix']+"_c1_grow_add.txt", params=PARAMS, gce=1, conf=1)
write_grow_file(filename=PARAMS['prefix']+"_grow_gibbs.txt", params=PARAMS, gce=2, conf=0, conf2=1)

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
MonteCarlo
RandomMT19937 seed {seed}
Configuration {vapor_config} particle_type0 {fstprt} cutoff {cutoff}
Configuration {liquid_config} particle_type0 {fstprt} cutoff {cutoff}
CopyFollowingLines for_num_configurations 2
    Potential Model LennardJones
    Potential Model LennardJones VisitModel VisitModelIntra intra_cut 3
    Potential VisitModel LongRangeCorrections
    RefPotential VisitModel DontVisitModel reference_index 0
EndCopy
ThermoParams beta {beta} chemical_potential 10
Metropolis
CopyNextLine replace0 configuration_index with0 0 replace1 tunable_param with1 30
TrialTranslate weight 0.5 tunable_param 1 configuration_index 1
CopyNextLine replace0 configuration_index with0 0 replace1 tunable_param with1 180
TrialParticlePivot weight 0.25 particle_type 0 tunable_param 0.4 pivot_site 0 configuration_index 1
CopyNextLine replace0 configuration_index with0 0 replace1 tunable_param with1 180
TrialParticlePivot weight 0.25 particle_type 0 tunable_param 0.4 pivot_site {last_site} configuration_index 1
TrialGrowFile grow_file {prefix}_c0_grow_canonical.txt
TrialGrowFile grow_file {prefix}_c1_grow_canonical.txt
CheckEnergy trials_per_update {tpc} decimal_places 4
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# gcmc initialization and nvt equilibration
CopyNextLine replace0 configuration_index with0 0 replace1 output_file with1 {prefix}{sim}_c0_fill.xyz
Movie trials_per_write {tpc} output_file {prefix}{sim}_c1_fill.xyz configuration_index 1
Log trials_per_write {tpc} output_file {prefix}{sim}_fill.csv include_bonds true
# decrease trials per due to infrequency of volume transfer attempts
Tune

# fill the first box
{init_vapor}

# fill the second box
{init_liquid}

Remove name0 Tune name1 Log name2 Movie name3 Movie

# gibbs equilibration cycles: equilibrate, estimate density, adjust, repeat
# start a very long run GibbsInitialize completes once targets are reached
Metropolis trials_per_cycle 1e9 cycles_to_complete 1e9
GibbsInitialize updates_density_equil {equil} updates_per_adjust {double_equil}
TrialGrowFile grow_file {prefix}_grow_gibbs.txt
TrialGibbsVolumeTransfer weight 0.006 tunable_param 3000 reference_index 0 print_num_accepted true
# a new tune is required when new Trials are introduced
Tune trials_per_tune 20
CheckEnergy trials_per_update {tpc} decimal_places 8
Log trials_per_write {tpc} output_file {prefix}{sim}_eq.csv
CopyNextLine replace0 configuration_index with0 0 replace1 output_file with1 {prefix}{sim}_c0_eq.xyz
Movie trials_per_write {tpc} output_file {prefix}{sim}_c1_eq.xyz configuration_index 1
ProfileCPU trials_per_write {tpc} output_file {prefix}{sim}_eq_profile.csv
# decrease trials per due to infrequency of volume transfer attempts
Run until complete
Remove name0 GibbsInitialize name1 Tune name2 Log name3 Movie name4 Movie name5 ProfileCPU

# gibbs ensemble production
Metropolis trials_per_cycle {tpc} cycles_to_complete {production_cycles}
Log trials_per_write {tpc} output_file {prefix}{sim}.csv
CopyFollowingLines for_num_configurations 2 replace_with_index [config]
    Movie   trials_per_write {tpc} output_file {prefix}{sim}_c[config].xyz
    Energy  trials_per_write {tpc} output_file {prefix}{sim}_c[config]_en.csv
    Density trials_per_write {tpc} output_file {prefix}{sim}_c[config]_dens.csv
EndCopy
GhostTrialVolume trials_per_update 1e3 trials_per_write {tpc} output_file {prefix}{sim}_pressure.csv
CPUTime trials_per_write {tpc} output_file {prefix}{sim}_cpu.csv
ProfileCPU trials_per_write {tpc} output_file {prefix}{sim}_profile.csv
Run until complete
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
    if len(diverged) > 0:
        print(diverged)
    assert len(diverged) == 0
    liquid_density = pd.read_csv(params['prefix']+"0_c1_dens.csv")
    liquid_density['average'] *= dens_conv
    liquid_density['block_stdev'] *= dens_conv
    liquid_density['diff'] = np.abs(liquid_density['average']-508)
    liquid_density['tol'] = np.sqrt(liquid_density['block_stdev']**2+(2**2))
    print(liquid_density)
    diverged = liquid_density[liquid_density['diff'] > z_factor*liquid_density['tol']]
    if len(diverged) > 0:
        print(diverged)
    assert len(diverged) == 0
    pres_conv = 1e33/na # convert from kJ/mol/A^3 to Pa (J/m^3)
    pressure = pd.read_csv(params['prefix']+"0_pressure.csv")
    pressure['average'] *= pres_conv
    pressure['block_stdev'] *= pres_conv
    pressure['diff'] = np.abs(pressure['average']-1.1976E+06)
    pressure['tol'] = np.sqrt(pressure['block_stdev']**2+(3.6212E+02)**2)
    print(pressure)

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
