import sys
import subprocess
import argparse
import random
import unittest

# define parameters of a pure component NVT MC hard sphere simulation
params = {
    "file_name": "post_process.xyz",
    "fstprt": "/feasst/forcefield/atom.fstprt",
    "num_minutes": 60,
    "num_nodes": 1, "procs_per_node": 1, "script": __file__}

# write fst script to run a single simulation
def mc_hs(params=params, file_name="launch.txt"):
    with open(file_name, "w") as myfile: myfile.write("""
MonteCarlo
Configuration particle_type0 {fstprt} xyz_file {file_name}
Potential Model HardSphere VisitModel VisitModelCell min_length 1
ThermoParams beta 1 chemical_potential 1
Metropolis
Scattering trials_per_update 1 trials_per_write 1 num_frequency 10 file_name hs_iq.csv
ReadConfigFromFile file_name {file_name}
Run until_criteria_complete true
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
        from pyfeasst import scattering
        iq=pd.read_csv('hs_iq.csv', comment="#")
        grp = iq.groupby('q', as_index=False)
        self.assertAlmostEqual(iq['i'][3810], 6.5772, delta=0.4)
        self.assertAlmostEqual(iq['i'][0]/iq['p0'][0]**2, 1, delta=0.075)

        import matplotlib.pyplot as plt
        plt.scatter(iq['q'], iq['i']/iq['p0']**2, label='sq_all')
        plt.plot(grp.mean()['q'], grp.mean()['i']/grp.mean()['p0']**2, label='sq_av')
        plt.legend()
        plt.show()

# run the simulation and, if complete, analyze.
def run():
    if args.task == 0:
        file_name = "hs_launch.txt"
        mc_hs(params, file_name=file_name)
        syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > hs_launch.log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../../../build/bin/rst hs_checkpoint.fst", shell=True, executable='/bin/bash')
    if syscode == 0:
        unittest.main(argv=[''], verbosity=2, exit=False)
    return syscode

if __name__ == "__main__":
    if args.run_type == 0:
        slurm_queue()
        subprocess.call("sbatch --array=0-10%1 slurm.txt | awk '{print $4}' >> launch_ids.txt", shell=True, executable='/bin/bash')
    elif args.run_type == 1:
        syscode = run()
        if syscode != 0:
            sys.exit(1)
    elif args.run_type == 2:
        unittest.main(argv=[''], verbosity=2, exit=False)
    else:
        assert False  # unrecognized run_type
