"""
This is a reproduction of the work described in https://doi.org/10.1016/j.xphs.2018.12.013
Group a mAb into domains and use the pdb file to compute the domain center of mass positions, bond lengths and angles.
Use Mayer-sampling simulations of individual domains to compute the excluded volume
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import feasstio
from pyfeasst import coarse_grain_pdb

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/particle/lj.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--pdb_file', type=str, default="../../../pyfeasst/tests/1igt.pdb",
                    help='pdb file that describes a mAb')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e5),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=int(1e0),
                    help='number of iterations for equilibration')
PARSER.add_argument('--production_iterations', type=int, default=int(1e1),
                    help='number of iterations for production')
PARSER.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=9, help='number of processors')
PARSER.add_argument('--prefix', type=str, default='cg', help='prefix for all output file names')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--seed', type=int, default=-1,
                    help='Random number generator seed. If -1, assign random seed to each sim.')
PARSER.add_argument('--max_restarts', type=int, default=0, help='Number of restarts in queue')
PARSER.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
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
PARAMS['script'] = __file__
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.99*PARAMS['hours_terminate'] - 0.0333 # terminate before queue
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']
def sim_node_dependent_params(params):
    """ Set parameters that depent upon the sim or node here. """
    if params['sim'] == 0: params['domain'] = 'fc'
    if params['sim'] == 1: params['domain'] = 'fab1'
    if params['sim'] == 2: params['domain'] = 'fab2'
    if params['sim'] == 3: params['domain'] = 'fv1'
    if params['sim'] == 4: params['domain'] = 'fv2'
    if params['sim'] == 5: params['domain'] = 'ch1_1'
    if params['sim'] == 6: params['domain'] = 'ch1_2'
    if params['sim'] == 7: params['domain'] = 'ch2'
    if params['sim'] == 8: params['domain'] = 'ch3'

# From table S2 of https://doi.org/10.1016/j.xphs.2018.12.013
# Heavy chains are B and D, while light chains are A and C, for fab1 and fab2, respectively.
# note that PDB coordinates are in Angstroms, while the manuscript is in nanometers
chains = {
          'hinge': {'B': range(236, 244), 'D': range(236, 244)},
          'fc': {'B': range(248, 475), 'D': range(248, 475)},
          'fab1': {'A': range(1, 215), 'B': range(1, 230)},
          'fab2': {'C': range(1, 215), 'D': range(1, 230)},
          'fv1': {'A': range(1, 109), 'B': range(1, 113)},
          'fv2': {'C': range(1, 109), 'D': range(1, 113)},
          'ch1_1': {'A': range(109, 215), 'B': range(113, 230)},
          'ch1_2': {'C': range(109, 215), 'D': range(113, 230)},
          'ch2': {'B': range(248, 361), 'D': range(248, 361)},
          'ch3': {'B': range(361, 475), 'D': range(361, 475)}}

# 4 bead (fab1, fab2, fc and hinge)
fc = coarse_grain_pdb.subset(pdb_file=ARGS.pdb_file, chains=chains['fc'])
r_com_fc = coarse_grain_pdb.center_of_mass(fc)/10  # divide all COM by 10 for Angstrom to nm
hinge = coarse_grain_pdb.subset(pdb_file=ARGS.pdb_file, chains=chains['hinge'])
r_com_hinge = coarse_grain_pdb.center_of_mass(hinge)/10
fab1 = coarse_grain_pdb.subset(pdb_file=ARGS.pdb_file, chains=chains['fab1'])
r_com_fab1 = coarse_grain_pdb.center_of_mass(fab1)/10
fab2 = coarse_grain_pdb.subset(pdb_file=ARGS.pdb_file, chains=chains['fab2'])
r_com_fab2 = coarse_grain_pdb.center_of_mass(fab2)/10

coarse_grain_pdb.pdb_to_fstprt(hinge, '1igt_hinge.fstprt')
coarse_grain_pdb.pdb_to_fstprt(fc, '1igt_fc.fstprt')
coarse_grain_pdb.pdb_to_fstprt(fab1, '1igt_fab1.fstprt')
coarse_grain_pdb.pdb_to_fstprt(fab2, '1igt_fab2.fstprt')

# 7 bead (fv[1,2], ch1_[1,2], ch2, ch3 and hinge)
fv1 = coarse_grain_pdb.subset(pdb_file=ARGS.pdb_file, chains=chains['fv1'])
r_com_fv1 = coarse_grain_pdb.center_of_mass(fv1)/10
fv2 = coarse_grain_pdb.subset(pdb_file=ARGS.pdb_file, chains=chains['fv2'])
r_com_fv2 = coarse_grain_pdb.center_of_mass(fv2)/10
ch1_1 = coarse_grain_pdb.subset(pdb_file=ARGS.pdb_file, chains=chains['ch1_1'])
r_com_ch1_1 = coarse_grain_pdb.center_of_mass(ch1_1)/10
ch1_2 = coarse_grain_pdb.subset(pdb_file=ARGS.pdb_file, chains=chains['ch1_2'])
r_com_ch1_2 = coarse_grain_pdb.center_of_mass(ch1_2)/10
ch2 = coarse_grain_pdb.subset(pdb_file=ARGS.pdb_file, chains=chains['ch2'])
r_com_ch2 = coarse_grain_pdb.center_of_mass(ch2)/10
ch3 = coarse_grain_pdb.subset(pdb_file=ARGS.pdb_file, chains=chains['ch3'])
r_com_ch3 = coarse_grain_pdb.center_of_mass(ch3)/10

coarse_grain_pdb.pdb_to_fstprt(fv1, '1igt_fv1.fstprt')
coarse_grain_pdb.pdb_to_fstprt(fv2, '1igt_fv2.fstprt')
coarse_grain_pdb.pdb_to_fstprt(ch1_1, '1igt_ch1_1.fstprt')
coarse_grain_pdb.pdb_to_fstprt(ch1_2, '1igt_ch1_2.fstprt')
coarse_grain_pdb.pdb_to_fstprt(ch2, '1igt_ch2.fstprt')
coarse_grain_pdb.pdb_to_fstprt(ch3, '1igt_ch3.fstprt')

# compute the distances and angles between the COM of pairs and triplets of domains
# compare with table S1 of https://doi.org/10.1016/j.xphs.2018.12.013
fc_hinge = r_com_fc - r_com_hinge
d_fc_hinge = np.sqrt(np.dot(fc_hinge, fc_hinge))
print('fc-hinge', d_fc_hinge, 'nm vs 4.24')
fab1_hinge = r_com_fab1 - r_com_hinge
d_fab1_hinge = np.sqrt(np.dot(fab1_hinge, fab1_hinge))
print('fab1-hinge', d_fab1_hinge, 'nm vs 5.85')
fab2_hinge = r_com_fab2 - r_com_hinge
d_fab2_hinge = np.sqrt(np.dot(fab2_hinge, fab2_hinge))
print('fab2-hinge', d_fab2_hinge, 'nm vs 4.87')
ch2_hinge = r_com_ch2 - r_com_hinge
d_ch2_hinge = np.sqrt(np.dot(ch2_hinge, ch2_hinge))
print('ch2-hinge', d_ch2_hinge, 'nm vs 2.98')
ch1_1_hinge = r_com_ch1_1 - r_com_hinge
d_ch1_1_hinge = np.sqrt(np.dot(ch1_1_hinge, ch1_1_hinge))
print('ch1_1-hinge', d_ch1_1_hinge, 'nm vs 4.15')
ch1_2_hinge = r_com_ch1_2 - r_com_hinge
d_ch1_2_hinge = np.sqrt(np.dot(ch1_2_hinge, ch1_2_hinge))
print('ch1_2-hinge', d_ch1_2_hinge, 'nm vs 3.42')
ch2_ch3 = r_com_ch2 - r_com_ch3
d_ch2_ch3 = np.sqrt(np.dot(ch2_ch3, ch2_ch3))
print('ch2-ch3', d_ch2_ch3, 'nm vs 2.83')
fv1_ch1_1 = r_com_fv1 - r_com_ch1_1
d_fv1_ch1_1 = np.sqrt(np.dot(fv1_ch1_1, fv1_ch1_1))
print('fv1-ch1_1', d_fv1_ch1_1, 'nm vs 3.48')
fv2_ch1_2 = r_com_fv2 - r_com_ch1_2
d_fv2_ch1_2 = np.sqrt(np.dot(fv2_ch1_2, fv2_ch1_2))
print('fv2-ch1_2', d_fv2_ch1_2, 'nm vs 3.15')
print('fc-h-fab1', np.arccos(np.dot(fc_hinge, fab1_hinge)/d_fc_hinge/d_fab1_hinge)*180/np.pi, 'degrees vs 107.17')
print('fc-h-fab2', np.arccos(np.dot(fc_hinge, fab2_hinge)/d_fc_hinge/d_fab2_hinge)*180/np.pi, 'degrees vs 108.98')
print('fab1-h-fab2', np.arccos(np.dot(fab1_hinge, fab2_hinge)/d_fab1_hinge/d_fab2_hinge)*180/np.pi, 'degrees vs 122.58')
print('ch2-h-ch1_1', np.arccos(np.dot(ch2_hinge, ch1_1_hinge)/d_ch2_hinge/d_ch1_1_hinge)*180/np.pi, 'degrees vs 122.31')
print('ch2-h-ch1_2', np.arccos(np.dot(ch2_hinge, ch1_2_hinge)/d_ch2_hinge/d_ch1_2_hinge)*180/np.pi, 'degrees vs 107.15')
print('ch1_1-h-ch1_2', np.arccos(np.dot(ch1_1_hinge, ch1_2_hinge)/d_ch1_1_hinge/d_ch1_2_hinge)*180/np.pi, 'degrees vs 98.30')

# compute the hinge radius of gyration
# compare with table 1 of https://doi.org/10.1016/j.xphs.2018.12.013
x_c = np.average(hinge['x_coord'])
y_c = np.average(hinge['y_coord'])
z_c = np.average(hinge['z_coord'])
rg2 = 0
for index, x in enumerate(hinge['x_coord']):
    dy = hinge['y_coord'].values[index] - y_c
    dz = hinge['z_coord'].values[index] - z_c
    rg2 += ((x-x_c)*(x-x_c) + dy*dy + dz*dz)
rg2 /= len(hinge['x_coord'])
print('2rg=sigma_hinge', 2*np.sqrt(rg2)/10, 'nm vs 1.52')

def write_feasst_script(params, file_name):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(file_name, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length 200 particle_type0 1igt_{domain}.fstprt \
    add_particles_of_type0 2 \
    group0 first first_particle_index 0 \
    group1 com com_site_type 5
Potential Model HardSphere VisitModel VisitModelCell min_length 3.9 energy_cutoff 1e100
RefPotential Model HardSphere sigma0 0 sigma1 0 sigma2 0 sigma3 0 sigma4 0 sigma5 30 cutoff0 0 cutoff1 0 cutoff2 0 cutoff3 0 cutoff4 0 cutoff5 30 group com
ThermoParams beta 1
MayerSampling num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
TrialTranslate new_only true reference_index 0 tunable_param 1 group first
TrialRotate new_only true reference_index 0 tunable_param 40
Checkpoint file_name {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}
set_variable trials_per 1e4

# tune trial parameters
CriteriaWriter trials_per_write trials_per file_name {prefix}_{domain}_b2_eq.txt
#Log trials_per_write trials_per file_name {prefix}_{domain}_eq.txt
#Movie trials_per_write trials_per file_name {prefix}_{domain}_eq.xyz
Tune
Run until_criteria_complete true
RemoveModify name Tune

# production
CriteriaWriter trials_per_write trials_per file_name {prefix}_{domain}_b2.txt
#Log trials_per_write trials_per file_name {prefix}_{domain}.txt
#Movie trials_per_write trials_per file_name {prefix}_{domain}.xyz
MayerSampling num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {production_iterations}
Run until_criteria_complete true
""".format(**params))

def post_process(params):
    def b2(file_name):
        file1 = open(file_name, 'r')
        lines = file1.readlines()
        file1.close()
        exec('iprm=' + lines[0], globals())
        return iprm
    b2hs_ref = 2*np.pi*3**3/3 # sigma=3 nanometer reference HS
    fc = b2(params['prefix']+'_fc_b2.txt')
    print('fc', fc['second_virial_ratio']*b2hs_ref, '+/-', fc['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 527.87 ± 1.91')
    fab1 = b2(params['prefix']+'_fab1_b2.txt')
    print('fab1', fab1['second_virial_ratio']*b2hs_ref, '+/-', fab1['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 443.20 ± 0.26')
    fab2 = b2(params['prefix']+'_fab2_b2.txt')
    print('fab2', fab2['second_virial_ratio']*b2hs_ref, '+/-', fab2['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 443.20 ± 0.26')
    fv1 = b2(params['prefix']+'_fv1_b2.txt')
    print('fv1', fv1['second_virial_ratio']*b2hs_ref, '+/-', fv1['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 208.13 ± 018')
    fv2 = b2(params['prefix']+'_fv2_b2.txt')
    print('fv2', fv2['second_virial_ratio']*b2hs_ref, '+/-', fv2['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 208.13 ± 018')
    ch1_1 = b2(params['prefix']+'_ch1_1_b2.txt')
    print('ch1_1', ch1_1['second_virial_ratio']*b2hs_ref, '+/-', ch1_1['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 179.09 ± 0.06')
    ch1_2 = b2(params['prefix']+'_ch1_2_b2.txt')
    print('ch1_2', ch1_2['second_virial_ratio']*b2hs_ref, '+/-', ch1_2['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 179.09 ± 0.06')
    ch2 = b2(params['prefix']+'_ch2_b2.txt')
    print('ch2', ch2['second_virial_ratio']*b2hs_ref, '+/-', ch2['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 316.83 ± 0.62')
    ch3 = b2(params['prefix']+'_ch3_b2.txt')
    print('ch3', ch3['second_virial_ratio']*b2hs_ref, '+/-', ch3['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 196.05 ± 0.14')

if __name__ == '__main__':
    feasstio.run_simulations(params=PARAMS,
                             sim_node_dependent_params=sim_node_dependent_params,
                             write_feasst_script=write_feasst_script,
                             post_process=post_process,
                             queue_function=feasstio.slurm_single_node,
                             args=ARGS)
