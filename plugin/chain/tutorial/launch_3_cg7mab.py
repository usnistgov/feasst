import sys
import subprocess
import argparse
import random
import unittest
import numpy as np
from multiprocessing import Pool

params = {
    "cubic_box_length": 90, "fstprt": "/feasst/plugin/chain/forcefield/cg7mab2.fstprt",
    "trials_per": 1e5, "hours_per_checkpoint": 1, "seed": random.randrange(int(1e9)), "num_hours": 0.8,
    "equilibration": 1e5, "production": 1e7, "num_nodes": 1, "procs_per_node": 1, "script": __file__}
params["num_minutes"] = round(params["num_hours"]*60)
params["num_hours_terminate"] = 0.95*params["num_hours"]*params["procs_per_node"]
nums=[30,3]
#nums=[605,530,454,378,303,227,151,76,30,15,6,3]

# write fst script to run a single simulation
def mc_cg7mab(params=params, file_name="launch.txt"):
    with open(file_name, "w") as myfile: myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt}
Potential Model HardSphere VisitModel VisitModelCell min_length 5.3283
ThermoParams beta 1 chemical_potential 1
Metropolis
TrialTranslate weight 0.5 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 0.5 tunable_param 0.2 tunable_target_acceptance 0.25 particle_type 0
Log trials_per_write {trials_per} file_name cg7mab_n{num_particles}.csv
Tune
CheckEnergy trials_per_update {trials_per} tolerance 1e-8
Checkpoint file_name cg7mab_n{num_particles}.fst num_hours_terminate {num_hours_terminate}

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Run until_num_particles {num_particles}
RemoveTrial name TrialAdd
Run num_trials {equilibration}

# nvt production
Movie trials_per_write {trials_per} file_name cg7mab_n{num_particles}.xyz
PairDistribution trials_per_update 1000 trials_per_write {trials_per} dr 0.025 file_name cg7mab_gr_n{num_particles}.csv print_intra true
Scattering trials_per_update 100 trials_per_write {trials_per} num_frequency 4 file_name cg7mab_iq_n{num_particles}.csv
Run num_trials {production}
""".format(**params))

# write slurm script
def slurm_queue():
    with open("slurm.txt", "w") as myfile: myfile.write("""#!/bin/bash
#SBATCH -n {procs_per_node} -N {num_nodes} -t {num_minutes}:00 -o hostname_%j.out -e hostname_%j.out
echo "Running {script} ID $SLURM_JOB_ID on $(hostname) at $(date) in $PWD"
cd $PWD
export OMP_NUM_THREADS={procs_per_node}
python {script} --run_type 1 --task $SLURM_ARRAY_TASK_ID
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
parser.add_argument('--task', type=int, default=0, help="input by slurm scheduler. If >0, restart from checkpoint.")
args = parser.parse_args()

# after the simulation is complete, perform some analysis
class TestFlatHistogramHS(unittest.TestCase):
    def test(self):
        import numpy as np
        import math
        import pandas as pd
        import matplotlib.pyplot as plt
        from pyfeasst import scattering
        #gr3=pd.read_csv('cg7mab_gr_n3.csv', comment="#")
        #gr30=pd.read_csv('cg7mab_gr_n30.csv', comment="#")
        iq3=pd.read_csv('cg7mab_iq_n3.csv', comment="#")
        iq30=pd.read_csv('cg7mab_iq_n30.csv', comment="#")
        grp3 = iq3.groupby('q', as_index=False)
        grp30 = iq30.groupby('q', as_index=False)
        plt.scatter(grp3.mean()['q'], grp30.mean()['i']/grp3.mean()['i'], label='direct I(10g/L)/I(1g/L)', marker='.')
        iq3rdfft = scattering.intensity(gr_file='cg7mab_gr_n3.csv', iq_file='cg7mab_iq_n3.csv', num_density=3/90**3, skip=10)
        iq30rdfft = scattering.intensity(gr_file='cg7mab_gr_n30.csv', iq_file='cg7mab_iq_n30.csv', num_density=30/90**3, skip=10)
        plt.scatter(iq3rdfft['q'], iq30rdfft['iq']/iq3rdfft['iq'], label='rdf ft I(10g/L)/I(1g/L)')
#        plt.ylabel('S', fontsize=16)
#        plt.xlabel('q(1/nm)', fontsize=16)
#        plt.legend()
#        plt.show()
        self.assertAlmostEqual(grp30.mean()['i'][0]/grp3.mean()['i'][0], 0.82777, delta=0.05)
        self.assertAlmostEqual(iq30rdfft['iq'][0]/iq3rdfft['iq'][0], 0.7784, delta=0.05)

# run the simulation and, if complete, analyze.
def run(nump):
    params["num_particles"] = nump
    if args.task == 0:
        file_name = "cg7mab_launch"+str(nump)+".txt"
        mc_cg7mab(params, file_name=file_name)
        syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > cg7mab_launch"+str(nump)+".log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../../../build/bin/rst cg7mab_n" + str(nump) + ".fst", shell=True, executable='/bin/bash')
    if syscode == 0:
        unittest.main(argv=[''], verbosity=2, exit=False)
    return syscode

if __name__ == "__main__":
    if args.run_type == 0:
        slurm_queue()
        subprocess.call("sbatch --array=0-10%1 slurm.txt | awk '{print $4}' >> launch_ids.txt", shell=True, executable='/bin/bash')
    elif args.run_type == 1:
        with Pool(len(nums)) as pool:
            codes = pool.starmap(run, zip(nums))
            if np.count_nonzero(codes) > 0:
                sys.exit(1)
    else:
        assert False  # unrecognized run_type
