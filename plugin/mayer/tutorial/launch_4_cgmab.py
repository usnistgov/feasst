# This is a reproduction of the work described in https://doi.org/10.1016/j.xphs.2018.12.013
# Group a mAb into domains and use the pdb file to compute the domain center of mass positions, bond lengths and angles.
# Use Mayer-sampling simulations of individual domains to compute the excluded volume
import random
import unittest
import argparse
import sys
import subprocess
import numpy as np
from multiprocessing import Pool
from pyfeasst import coarse_grain_pdb

pdb_file = "../../../pyfeasst/tests/1igt.pdb"

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
fc = coarse_grain_pdb.subset(pdb_file=pdb_file, chains=chains['fc'])
r_com_fc = coarse_grain_pdb.center_of_mass(fc)/10  # divide all COM by 10 for Angstrom to nm
hinge = coarse_grain_pdb.subset(pdb_file=pdb_file, chains=chains['hinge'])
r_com_hinge = coarse_grain_pdb.center_of_mass(hinge)/10
fab1 = coarse_grain_pdb.subset(pdb_file=pdb_file, chains=chains['fab1'])
r_com_fab1 = coarse_grain_pdb.center_of_mass(fab1)/10
fab2 = coarse_grain_pdb.subset(pdb_file=pdb_file, chains=chains['fab2'])
r_com_fab2 = coarse_grain_pdb.center_of_mass(fab2)/10

coarse_grain_pdb.pdb_to_fstprt(hinge, '1igt_hinge.fstprt')
coarse_grain_pdb.pdb_to_fstprt(fc, '1igt_fc.fstprt')
coarse_grain_pdb.pdb_to_fstprt(fab1, '1igt_fab1.fstprt')
coarse_grain_pdb.pdb_to_fstprt(fab2, '1igt_fab2.fstprt')

# 7 bead (fv[1,2], ch1_[1,2], ch2, ch3 and hinge)
fv1 = coarse_grain_pdb.subset(pdb_file=pdb_file, chains=chains['fv1'])
r_com_fv1 = coarse_grain_pdb.center_of_mass(fv1)/10
fv2 = coarse_grain_pdb.subset(pdb_file=pdb_file, chains=chains['fv2'])
r_com_fv2 = coarse_grain_pdb.center_of_mass(fv2)/10
ch1_1 = coarse_grain_pdb.subset(pdb_file=pdb_file, chains=chains['ch1_1'])
r_com_ch1_1 = coarse_grain_pdb.center_of_mass(ch1_1)/10
ch1_2 = coarse_grain_pdb.subset(pdb_file=pdb_file, chains=chains['ch1_2'])
r_com_ch1_2 = coarse_grain_pdb.center_of_mass(ch1_2)/10
ch2 = coarse_grain_pdb.subset(pdb_file=pdb_file, chains=chains['ch2'])
r_com_ch2 = coarse_grain_pdb.center_of_mass(ch2)/10
ch3 = coarse_grain_pdb.subset(pdb_file=pdb_file, chains=chains['ch3'])
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

# Define a Mayer-sampling simulation using FEASST
def mc(file_name):
    with open(file_name, 'w') as file: file.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_box_length 200 particle_type0 1igt_{domain}.fstprt \
    add_particles_of_type0 2 \
    group0 first first_particle_index 0 \
    group1 com com_site_type 5
Potential Model HardSphere VisitModel VisitModelCell min_length 3.9 energy_cutoff 1e100
RefPotential Model HardSphere sigma0 0 sigma1 0 sigma2 0 sigma3 0 sigma4 0 sigma5 30 cutoff0 0 cutoff1 0 cutoff2 0 cutoff3 0 cutoff4 0 cutoff5 30 group com
ThermoParams beta 1
MayerSampling
TrialTranslate new_only true reference_index 0 tunable_param 1 group first
TrialRotate new_only true reference_index 0 tunable_param 40
set_variable trials_per 1e4

# tune trial parameters
CriteriaWriter trials_per trials_per file_name cg_b2_eq_{domain}.txt
#Log trials_per trials_per file_name cg_eq_{domain}.txt
#Movie trials_per trials_per file_name cg_eq_{domain}.xyz
Tune
Run num_trials 1e5
RemoveModify name Tune

# production
CriteriaWriter trials_per trials_per file_name cg_b2_{domain}.txt
#Log trials_per trials_per file_name cg_{domain}.txt
#Movie trials_per trials_per file_name cg_{domain}.xyz
MayerSampling
Run num_trials 1e6
""".format(**params))

# write slurm script to fill an HPC node with simulations
def slurm_queue():
    with open("slurm.txt", "w") as myfile: myfile.write("""#!/bin/bash
#SBATCH -n {procs_per_node} -N 1 -t 1440:00 -o hostname_%j.out -e hostname_%j.out
echo "Running ID $SLURM_JOB_ID on $(hostname) at $(date) in $PWD"
cd $PWD
python {script} --run_type 1
if [ $? == 0 ]; then
  echo "Job is done"
  scancel $SLURM_ARRAY_JOB_ID
else
  echo "Job is terminating, to be restarted again"
fi
echo "Time is $(date)"
""".format(**params))

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--run_type', '-r', type=int, default=0, help="0: submit batch to scheduler, 1: run batch on host")
args = parser.parse_args()

params={'procs_per_node': 9, 'script': __file__}

# run a single simulation as part of the batch to fill a node
def run(proc):
    if proc == 0: params['domain'] = 'fc'
    if proc == 1: params['domain'] = 'fab1'
    if proc == 2: params['domain'] = 'fab2'
    if proc == 3: params['domain'] = 'fv1'
    if proc == 4: params['domain'] = 'fv2'
    if proc == 5: params['domain'] = 'ch1_1'
    if proc == 6: params['domain'] = 'ch1_2'
    if proc == 7: params['domain'] = 'ch2'
    if proc == 8: params['domain'] = 'ch3'
    params["seed"] = random.randrange(int(1e9))
    file_name = "launch_run"+params['domain']+".txt"
    mc(file_name)
    syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > launch_run"+params['domain']+".log", shell=True, executable='/bin/bash')
    if syscode == 0:
        unittest.main(argv=[''], verbosity=2, exit=False)
    return syscode

# after the simulation is complete, perform some analysis
class TestCGmAb(unittest.TestCase):
    def test(self):
        def b2(file_name):
            file1 = open(file_name, 'r')
            lines = file1.readlines()
            file1.close()
            exec('iprm=' + lines[0], globals())
            return iprm
        b2hs_ref = 2*np.pi*3**3/3 # sigma=3 nanometer reference HS
        fc = b2('cg_b2_fc.txt')
        print('fc', fc['second_virial_ratio']*b2hs_ref, '+/-', fc['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 527.87 ± 1.91')
        fab1 = b2('cg_b2_fab1.txt')
        print('fab1', fab1['second_virial_ratio']*b2hs_ref, '+/-', fab1['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 443.20 ± 0.26')
        fab2 = b2('cg_b2_fab2.txt')
        print('fab2', fab2['second_virial_ratio']*b2hs_ref, '+/-', fab2['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 443.20 ± 0.26')
        fv1 = b2('cg_b2_fv1.txt')
        print('fv1', fv1['second_virial_ratio']*b2hs_ref, '+/-', fv1['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 208.13 ± 018')
        fv2 = b2('cg_b2_fv2.txt')
        print('fv2', fv2['second_virial_ratio']*b2hs_ref, '+/-', fv2['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 208.13 ± 018')
        ch1_1 = b2('cg_b2_ch1_1.txt')
        print('ch1_1', ch1_1['second_virial_ratio']*b2hs_ref, '+/-', ch1_1['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 179.09 ± 0.06')
        ch1_2 = b2('cg_b2_ch1_2.txt')
        print('ch1_2', ch1_2['second_virial_ratio']*b2hs_ref, '+/-', ch1_2['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 179.09 ± 0.06')
        ch2 = b2('cg_b2_ch2.txt')
        print('ch2', ch2['second_virial_ratio']*b2hs_ref, '+/-', ch2['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 316.83 ± 0.62')
        ch3 = b2('cg_b2_ch3.txt')
        print('ch3', ch3['second_virial_ratio']*b2hs_ref, '+/-', ch3['second_virial_ratio_block_stdev']*b2hs_ref, 'nm^3 vs 196.05 ± 0.14')

if __name__ == "__main__":
    if args.run_type == 0:
        slurm_queue()
        subprocess.call("sbatch slurm.txt | awk '{print $4}' >> launch_ids.txt", shell=True, executable='/bin/bash')
    elif args.run_type == 1:
        with Pool(params["procs_per_node"]) as pool:
            codes = pool.starmap(run, zip(range(0, params["procs_per_node"])))
            if np.count_nonzero(codes) > 0:
                sys.exit(1)
    else:
        assert False  # unrecognized run_type
