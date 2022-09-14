import sys
import subprocess
import numpy as np
import argparse
from multiprocessing import Pool
import random
import unittest

# define parameters of a pure component NVT MC Lennard-Jones simulation
params = {
    "num_particles": 500, "density": 0.001, "trials_per": 1e5,
    "beta": 1./0.9, "fstprt": "/feasst/forcefield/lj.fstprt",
    "equilibration": 1e7, "production": 1e7,
    "seed": random.randrange(int(1e9)), "num_hours": 1, "script": __file__}
params["box_length"] = (params["num_particles"]/params["density"])**(1./3.)
params["num_minutes"] = round(params["num_hours"]*60)
params["num_hours_terminate"] = 0.95*params["num_hours"]

# write fst script to run a single simulation
def mc_lj(params=params, file_name="launch.txt"):
    with open(file_name, "w") as myfile: myfile.write("""
MonteCarlo
Checkpoint file_name checkpoint.fst
RandomMT19937 seed time
Configuration cubic_box_length {box_length} particle_type0 {fstprt}
Potential Model LennardJones VisitModel VisitModelCell min_length 3
Potential VisitModel LongRangeCorrections
ThermoParams beta 0.1 chemical_potential 10
Metropolis
TrialTranslate tunable_param 2. tunable_target_acceptance 0.2
TrialAdd particle_type 0
Run until_num_particles {num_particles}
RemoveTrial name TrialAdd
ThermoParams beta {beta}
Tune
CheckEnergy trials_per {trials_per} tolerance 1e-8

# equilibrate
Run num_trials {equilibration}
RemoveModify name Tune

# production analysis and output
Log trials_per {trials_per} file_name lj.txt
Energy trials_per_write {trials_per} file_name en.txt
Run num_trials {production}
""".format(**params))

# write slurm script to fill nodes with simulations
def slurm_queue():
    with open("slurm.txt", "w") as myfile: myfile.write("""#!/bin/bash
#SBATCH -n {procs_per_node} -N {num_nodes} -t {num_minutes}:00 -o hostname_%j.out -e hostname_%j.out
echo "Running ID $SLURM_JOB_ID on $(hostname) at $(date) in $PWD"
cd $PWD
python {script} --run_type 1 --task $SLURM_ARRAY_TASK_ID
if [ $? == 0 ]; then
  echo "Job is done"
  scancel $SLURM_ARRAY_JOB_ID
else
  echo "Job is terminating, to be restarted again"
fi
echo "Time is $(date)"
""".format(**params))

# add additional arguments for multi-core simulations
params.update({"sim": 0, "num_nodes": 1, "procs_per_node": 1})
params["num_sims"] = params["num_nodes"]*params["procs_per_node"]

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--run_type', '-r', type=int, default=0, help="0: submit batch to scheduler, 1: run batch on host")
parser.add_argument('--task', type=int, default=0, help="input by slurm scheduler. If >0, restart from checkpoint.")
args = parser.parse_args()

# after the simulation is complete, perform some tests
class TestMonteCarloLJ(unittest.TestCase):
    def test(self):
        # test the average energy against the NIST SRSW
        import pandas as pd
        df = pd.read_csv('en.txt')
        stdev = (df['block_stdev'][0]**2 + (1.89E-05)**2)**(1./2.)
        self.assertAlmostEqual(-9.9165E-03*params["num_particles"], df['average'][0], delta=2.576*stdev)

# run a single simulation as part of the batch to fill a node
def run(sim):
    if args.task == 0:
        params["sim"] = sim
        params["seed"] = random.randrange(int(1e9))
        file_name = "launch_run"+str(sim)+".txt"
        mc_lj(params, file_name=file_name)
        syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > launch_run"+str(sim)+".log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../../../build/bin/rst checkpoint" + str(sim) + ".fst", shell=True, executable='/bin/bash')
    if syscode == 0:
        unittest.main(argv=[''], verbosity=2, exit=False)
    return syscode

if __name__ == "__main__":
    if args.run_type == 0:
        slurm_queue()
        subprocess.call("sbatch --array=0-10%1 slurm.txt | awk '{print $4}' >> launch_ids.txt", shell=True, executable='/bin/bash')
    elif args.run_type == 1:
        with Pool(params["num_sims"]) as pool:
            codes = pool.starmap(run, zip(range(0, params["num_sims"])))
            if np.count_nonzero(codes) > 0:
                sys.exit(1)
    else:
        assert False  # unrecognized run_type
