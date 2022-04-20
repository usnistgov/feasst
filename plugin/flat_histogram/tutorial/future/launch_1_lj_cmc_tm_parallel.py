import sys
import subprocess
import numpy as np
import argparse
from multiprocessing import Pool
import random

# define parameters of a pure component NVT MC Lennard-Jones simulation
params = {
    "cubic_box_length": 8, "fstprt": "/feasst/forcefield/lj.fstprt", "beta": 1/1.2,
    "max_particles": 370, "min_particles": 0, "min_sweeps": 1e5, "mu": -1.9,
    "trials_per": 1e5, "seed": random.randrange(1e9), "num_hours": 1,
    "equilibration": 1e6, "num_nodes": 1, "procs_per_node": 4}
params["num_minutes"] = round(params["num_hours"]*60)
params["num_hours"] = params["num_hours"]*params["procs_per_node"]
params["num_hours_terminate"] = 0.95*params["num_hours"]*params["procs_per_node"]

# write fst script to run a single simulation
def mc_lj(params=params, file_name="launch.txt"):
    with open(file_name, "w") as myfile: myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice min_window_size 8 hours_per 0.005 ln_prob_file lnpi.txt
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 1.5
Checkpoint file_name checkpoint.fst num_hours 0.1 num_hours_terminate {num_hours_terminate}

# begin description of each MC clone
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt}
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential {mu}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.2
Log trials_per {trials_per} file_name lj[sim_index].txt
Tune trials_per {trials_per}
CheckEnergy trials_per {trials_per} tolerance 1e-8
Checkpoint file_name checkpoint[sim_index].fst num_hours 0.1

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Run until_num_particles [soft_macro_min]
RemoveTrial name TrialAdd
Run num_trials {equilibration}

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
              Bias WLTM min_sweeps {min_sweeps} new_sweep 1 min_flatness 25 collect_flatness 20
#Bias TransitionMatrix min_sweeps {min_sweeps} new_sweep 1
TrialTransfer weight 2 particle_type 0
Movie trials_per {trials_per} file_name lj[sim_index].xyz
Energy trials_per_write {trials_per} file_name en[sim_index].txt multistate true
CriteriaUpdater trials_per {trials_per}
CriteriaWriter trials_per {trials_per} file_name crit[sim_index].txt
#Run until_criteria_complete true
""".format(**params))

# write slurm script to fill nodes with simulations
def slurm_queue():
    with open("slurm.txt", "w") as myfile: myfile.write("""#!/bin/bash
#SBATCH -n {procs_per_node} -N {num_nodes} -t {num_minutes}:00 -o hostname_%j.out -e hostname_%j.out
echo "Running ID $SLURM_JOB_ID on $(hostname) at $(date) in $PWD"
cd $PWD
python launch.py --run_type 1 --task $SLURM_ARRAY_TASK_ID
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

# run a single simulation as part of the batch to fill a node
def run():
    if args.task == 0:
        file_name = "launch_run.txt"
        mc_lj(params, file_name=file_name)
        syscode = subprocess.call("~/feasst/build/bin/fst < " + file_name + " > launch_run.log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("~/feasst/build/bin/rst checkpoint.fst", shell=True, executable='/bin/bash')
    return syscode

if __name__ == "__main__":
    if args.run_type == 0:
        slurm_queue()
        subprocess.call("sbatch --array=0-10%1 slurm.txt | awk '{print $4}' >> launch_ids.txt", shell=True, executable='/bin/bash')
    elif args.run_type == 1:
        syscode = run()
        if syscode > 0:
            sys.exit(1)
    else:
        assert(False) # unrecognized run_type
