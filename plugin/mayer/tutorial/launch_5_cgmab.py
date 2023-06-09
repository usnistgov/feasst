import random
import unittest
import argparse
import sys
import subprocess
import numpy as np
from multiprocessing import Pool
from pyfeasst import physical_constants
from pyfeasst import coarse_grain_pdb

params = {
    "seed": random.randrange(int(1e9)),
    "reference_sigma": 4.5,
    "trials_per": 1e4,
    'procs_per_node': 1, 'script': __file__}

def mc(file_name):
    with open(file_name, 'w') as file: file.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_box_length 500 particle_type0 /feasst/plugin/chain/forcefield/cg7mab2.fstprt \
    add_particles_of_type0 2 \
    group0 first first_particle_index 0
Potential Model HardSphere
#RefPotential Model HardSphere sigma0 {reference_sigma} sigma1 0 sigma2 0 sigma3 0 sigma4 0 sigma5 0 sigma6 0 \
#                       cutoff0 {reference_sigma} cutoff1 0 cutoff2 0 cutoff3 0 cutoff4 0 cutoff5 0 cutoff6 0
RefPotential Model HardSphere sigma 0 sigma0 {reference_sigma} cutoff 0 cutoff0 {reference_sigma}
ThermoParams beta 1
MayerSampling
TrialTranslate new_only true reference_index 0 tunable_param 1 group first
TrialRotate new_only true reference_index 0 tunable_param 40

# tune trial parameters
CriteriaWriter trials_per_write {trials_per} file_name cg_b2_eq.txt
Log trials_per_write {trials_per} file_name cg_eq.txt
Movie trials_per_write {trials_per} file_name cg_eq.xyz
Tune
Run num_trials 1e5
RemoveModify name Tune
RemoveAnalyze name CriteriaWriter

# production
CriteriaWriter trials_per_write {trials_per} file_name cg_b2.txt
Log trials_per_write {trials_per} file_name cg.txt
Movie trials_per_write {trials_per} file_name cg.xyz
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

# run a single simulation as part of the batch to fill a node
def run(proc):
    file_name = "launch_run.txt"
    mc(file_name)
    syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > launch_run.log", shell=True, executable='/bin/bash')
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
        b2hs_ref = 2*np.pi*params['reference_sigma']**3/3 # reference HS in nm^3
        b2 = b2('cg_b2.txt')
        print('b2', b2['second_virial_ratio']*b2hs_ref, '+/-', b2['second_virial_ratio_block_stdev']*b2hs_ref)
        lpm_fac = 1e-24*physical_constants.AvogadroConstant().value()
        print('b2', b2['second_virial_ratio']*b2hs_ref*lpm_fac,
             '+/-', b2['second_virial_ratio_block_stdev']*b2hs_ref*lpm_fac,
             'L/mol')
        mlmg_fac = 1e-24*physical_constants.AvogadroConstant().value()/150000
        print('b2', b2['second_virial_ratio']*b2hs_ref*mlmg_fac,
             '+/-', b2['second_virial_ratio_block_stdev']*b2hs_ref*mlmg_fac,
             'mL/mg')
        self.assertAlmostEqual(b2['second_virial_ratio']*b2hs_ref*mlmg_fac, 0.0113, delta=0.001)

if __name__ == "__main__":
    if args.run_type == 0:
        slurm_queue()
        subprocess.call("sbatch slurm.txt | awk '{print $4}' >> launch_ids.txt", shell=True, executable='/bin/bash')
    elif args.run_type == 1:
        with Pool(params["procs_per_node"]) as pool:
            codes = pool.starmap(run, zip(range(0, params["procs_per_node"])))
            if np.count_nonzero(codes) > 0:
                sys.exit(1)
    elif args.run_type == 2:
        unittest.main(argv=[''], verbosity=2, exit=False)
    else:
        assert False  # unrecognized run_type
