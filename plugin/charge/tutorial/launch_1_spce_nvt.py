# In this example, we reproduce the average energy reported in https://doi.org/10.1063/1.476834

import sys
import subprocess
import numpy as np
import argparse
from multiprocessing import Pool
import random
import unittest

# define parameters of a pure component NVT MC SPC/E simulation
params = {
    "num_particles": 512, "cubic_box_length": 24.8586887, "trials_per": 1e5,
    "temperature": 298, "fstprt": "/feasst/forcefield/spce.fstprt",
    "equilibration": 1e6, "production": 1e6,
    "seed": random.randrange(1e9), "num_hours": 1, "script": __file__}
params["cutoff"] = params["cubic_box_length"]/2.
params["alpha"] = 5.6/params["cubic_box_length"]
R = 1.3806488E-23*6.02214129E+23 # J/mol/K
params["beta"] = 1./(params["temperature"]*R/1e3) # mol/kJ
params["num_minutes"] = round(params["num_hours"]*60)
params["num_hours_terminate"] = 0.95*params["num_hours"]

# write fst script to run a single simulation
def mc_lj(params=params, file_name="launch.txt"):
    with open(file_name, "w") as myfile: myfile.write("""
MonteCarlo
Checkpoint file_name checkpoint.fst
RandomMT19937 seed time
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt} physical_constants CODATA2010 cutoff {cutoff}
Potential VisitModel Ewald alpha {alpha} kmax_squared 38
Potential Model ModelTwoBodyFactory model0 LennardJones model1 ChargeScreened VisitModel VisitModelCutoffOuter table_size 1e6
Potential Model ChargeScreenedIntra VisitModel VisitModelBond
Potential Model ChargeSelf
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential 1
Metropolis
TrialTranslate weight 0.5 tunable_param 0.275 tunable_target_acceptance 0.25
TrialParticlePivot weight 0.5 particle_type 0 tunable_param 0.5 tunable_target_acceptance 0.25
TrialAdd particle_type 0
Run until_num_particles {num_particles}
RemoveTrial name TrialAdd
ThermoParams beta {beta}
Log trials_per {trials_per} file_name spce.txt
Tune
CheckEnergy trials_per {trials_per} tolerance 1e-8

# equilibrate
Run num_trials {equilibration}
RemoveModify name Tune

# production analysis and output
Energy trials_per_write {trials_per} file_name spce_en.txt
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
        # test the average energy against https://doi.org/10.1063/1.476834
        # note that there are slight systematic differences in the energy due to different Ewald cutoffs, etc.
        import pandas as pd
        df = pd.read_csv('spce_en.txt')
        stdev = (df['block_stdev'][0]**2 + (0.02*params["num_particles"])**2)**(1./2.)
        self.assertAlmostEqual(-46.82*params["num_particles"], df['average'][0], delta=20*stdev)

# run a single simulation as part of the batch to fill a node
def run(sim):
    if args.task == 0:
        params["sim"] = sim
        params["seed"] = random.randrange(1e9)
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
