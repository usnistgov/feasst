import sys
import subprocess
import math
import argparse
import random
import unittest

# define parameters of a pure component MC Kern-Frenkel patchy particle fluid
# see https://doi.org/10.1063/1.1569473
params = {
    "cubic_box_length": 8, "fstprt": "/feasst/plugin/patch/forcefield/two_patch_linear.fstprt",
    "max_particles": 370, "min_particles": 0, "min_sweeps": 1e3, "mu": -1.5, "beta": 1/0.7, "chi": 0.7, "cutoff": 1.5,
    "trials_per": 1e5, "hours_per_adjust": 0.01, "hours_per_checkpoint": 1, "seed": random.randrange(1e9), "num_hours": 5*24,
    "equilibration": 1e5, "num_nodes": 1, "procs_per_node": 32, "script": __file__}
params["patch_angle"] = 2*math.asin(math.sqrt(params['chi']/2))*180/math.pi
params["num_minutes"] = round(params["num_hours"]*60)
params["hours_per_adjust"] = params["hours_per_adjust"]*params["procs_per_node"]
params["hours_per_checkpoint"] = params["hours_per_checkpoint"]*params["procs_per_node"]
params["num_hours_terminate"] = 0.95*params["num_hours"]*params["procs_per_node"]

# write fst script to run a single simulation
def mc_kf(params=params, file_name="launch.txt"):
    with open(file_name, "w") as myfile: myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_per_adjust} ln_prob_file kf_lnpi.txt bounds_file kf_bounds.txt num_adjust_per_write 10
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 2.5 min_size 2
Checkpoint file_name kf_checkpoint.fst num_hours {hours_per_checkpoint} num_hours_terminate {num_hours_terminate}

# begin description of each MC clone
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt} \
  patch_angle1 {patch_angle} group0 centers centers_site_type0 0 \
  cutoff {cutoff}
Potential Model HardSphere VisitModel VisitModelCell min_length 1 cell_group centers group centers
Potential Model SquareWell VisitModel VisitModelCell min_length {cutoff} cell_group centers \
  VisitModelInner VisitModelInnerPatch group centers
ThermoParams beta {beta} chemical_potential {mu}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
TrialRotate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
Log trials_per {trials_per} file_name kf[sim_index].txt
Tune
CheckEnergy trials_per {trials_per} tolerance 1e-8

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Run until_num_particles [soft_macro_min]
RemoveTrial name TrialAdd
Run num_trials {equilibration}
RemoveModify name Tune

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias TransitionMatrix min_sweeps {min_sweeps} new_sweep 1
TrialTransfer weight 2 particle_type 0
RemoveAnalyze name Log
Log trials_per {trials_per} file_name kf[sim_index].txt
Tune trials_per_write {trials_per} file_name kf_tune[sim_index].txt multistate true
MoviePatch trials_per {trials_per} file_name kf[sim_index].xyz
Energy trials_per_write {trials_per} file_name kf_en[sim_index].txt multistate true
CriteriaUpdater trials_per {trials_per}
CriteriaWriter trials_per {trials_per} file_name kf_crit[sim_index].txt
""".format(**params))

# write slurm script to fill nodes with simulations
def slurm_queue():
    with open("slurm.txt", "w") as myfile: myfile.write("""#!/bin/bash
#SBATCH -n {procs_per_node} -N {num_nodes} -t {num_minutes}:00 -o hostname_%j.out -e hostname_%j.out
echo "Running ID $SLURM_JOB_ID on $(hostname) at $(date) in $PWD"
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

# after the simulation is complete, perform some tests
class TestFlatHistogramKF(unittest.TestCase):
    def test(self):
        import numpy as np
        import pandas as pd
        lnpi=pd.read_csv('kf_lnpi.txt')
        lnpi=lnpi[:6] # cut down to three rows
        lnpi['ln_prob'] -= np.log(sum(np.exp(lnpi['ln_prob'])))  # renormalize
        lnpi['ln_prob_prev'] = [-15.9976474469475, -11.9104563420586,  -8.48324267323538, -5.42988602574393, -2.64984051640555, -0.07824246342703]
        lnpi['ln_prob_prev_stdev'] = [0.035, 0.03, 0.025, 0.02, 0.015, 0.01]
        diverged=lnpi[lnpi.ln_prob-lnpi.ln_prob_prev > 6*lnpi.ln_prob_prev_stdev]
        print(diverged)
        self.assertTrue(len(diverged) == 0)
        en=pd.read_csv('kf_en00.txt')
        en=en[:6]
        en['prev'] = [0, 0, -0.038758392176564, -0.116517384264731, -0.232665619265520, -0.387804181572135]
        diverged=en[en.average-en.prev > 10*en.block_stdev]
        print(diverged)
        self.assertTrue(len(diverged) == 0)

# run the simulation and, if complete, analyze.
def run():
    if args.task == 0:
        file_name = "kf_launch.txt"
        mc_kf(params, file_name=file_name)
        syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > kf_launch.log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../../../build/bin/rst kf_checkpoint.fst", shell=True, executable='/bin/bash')
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
    else:
        assert(False) # unrecognized run_type
