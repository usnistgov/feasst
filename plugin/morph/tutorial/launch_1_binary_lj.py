"""
Simulate Mixture I. of https://doi.org/10.1063/1.1844372
"""

import sys
import subprocess
import argparse
import random
import unittest
import pathlib
from pyfeasst import physical_constants

params = {
    "cubic_box_length": 7,
    "min_particles": 0,
    "beta": 0.8, "max_particles": 20, "mu0": -5.4, 'mu1': -5.5,
    "trials_per": 1e6, "hours_per_adjust": 0.01, "hours_per_checkpoint": 1, "seed": random.randrange(int(1e9)), "num_hours": 5*24,
    "equilibration": 1e6, "num_nodes": 1, "procs_per_node": 1, "script": __file__}
params["num_minutes"] = round(params["num_hours"]*60)
params["hours_per_adjust"] = params["hours_per_adjust"]*params["procs_per_node"]
params["hours_per_checkpoint"] = params["hours_per_checkpoint"]*params["procs_per_node"]
params["num_hours_terminate"] = 0.95*params["num_hours"]*params["procs_per_node"]
params["mu_init"] = -1
params["min_sweeps"]=1000
params["window_alpha"]=1.1
params["min_window_size"]=5

def write_fstprt():
    with open('lja.fstprt', 'w') as myfile: myfile.write("""
""".format(**params))


def mc_binary_lj(params, file_name):
    with open(file_name, "w") as myfile: myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_per_adjust} ln_prob_file binary_lj_lnpin{node}.txt bounds_file binary_lj_boundsn{node}.txt num_adjust_per_write 10 min_window_size {min_window_size}
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha {window_alpha} min_size {min_window_size}
Checkpoint file_name binary_lj_checkpointn{node}.fst num_hours {hours_per_checkpoint} num_hours_terminate {num_hours_terminate}

# begin description of each MC clone
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 /feasst/forcefield/atom.fstprt particle_type1 /feasst/forcefield/lj.fstprt sigma0 1.0 epsilon0 1.0 cutoff0 3.0 sigma1 1.064 epsilon1 1.37 cutoff1 3.0 sigma0_1 1.034 epsilon0_1 1.152 cutoff0_1 3.0
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential0 {mu_init} chemical_potential1 {mu_init}
Metropolis
TrialTranslate weight 0.5 tunable_param 0.2 tunable_target_acceptance 0.25
Log trials_per_write {trials_per} file_name binary_ljn{node}s[sim_index].txt
Tune
CheckEnergy trials_per_update {trials_per} tolerance 1e-4

# gcmc initialization and nvt equilibration
TrialAdd particle_type 1
Run until_num_particles {max_particles}
RemoveTrial name TrialAdd
ThermoParams beta {beta} chemical_potential0 {mu0} chemical_potential1 {mu1}
Metropolis
Run num_trials {equilibration}
RemoveModify name Tune

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles particle_type 0 width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} new_sweep 1 min_flatness 22 collect_flatness 20 min_collect_sweeps 18
TrialMorph particle_type0 0 particle_type_morph0 1
TrialMorph particle_type0 1 particle_type_morph0 0
RemoveAnalyze name Log
Log trials_per_write {trials_per} file_name binary_ljn{node}s[sim_index].txt
Movie trials_per_write {trials_per} file_name binary_ljn{node}s[sim_index].xyz
Tune trials_per_write {trials_per} file_name binary_lj_tunen{node}s[sim_index].txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per} file_name binary_lj_enn{node}s[sim_index].txt multistate true start_after_iteration 1
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {trials_per} file_name binary_lj_critn{node}s[sim_index].txt
""".format(**params))

# write slurm script
def slurm_queue(node):
    params["node"]=node
    with open("slurm"+str(node)+".txt", "w") as myfile: myfile.write("""#!/bin/bash
#SBATCH -n {procs_per_node} -N 1 -t {num_minutes}:00 -o hostname_%j.out -e hostname_%j.out
echo "Running {script} ID $SLURM_JOB_ID on $(hostname) at $(date) in $PWD"
cd $PWD
export OMP_NUM_THREADS={procs_per_node}
python {script} --run_type 1 --task $SLURM_ARRAY_TASK_ID --node {node}
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
parser.add_argument('--node', type=int, default=0, help="break the job into multiple nodes.")
args = parser.parse_args()
params['node']=args.node

# after the simulation is complete, perform some tests
class TestFlatHistogramBinaryLJ(unittest.TestCase):
    def test(self):
        import numpy as np
        import pandas as pd
        lnpi=pd.read_csv('binary_lj_lnpin0.txt')
        self.assertAlmostEqual(9.327631384282558, (np.exp(lnpi["ln_prob"]) * lnpi["state"]).sum(), delta=0.5)

# run the simulation and, if complete, analyze.
def run():
    if args.task == 0:
        file_name = "binary_lj_launch"+str(params["node"])+".txt"
        mc_binary_lj(params=params, file_name=file_name)
        syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > binary_lj_launch"+str(params['node'])+".log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../../../build/bin/rst binary_lj_checkpointn"+str(params['node'])+".fst", shell=True, executable='/bin/bash')
    if syscode == 0 and params['node'] == 1:
        unittest.main(argv=[''], verbosity=2, exit=False)
    return syscode

if __name__ == "__main__":
    if args.run_type == 0:
        for node in range(params['num_nodes']):
            slurm_queue(node)
            subprocess.call("sbatch --array=0-2%1 slurm"+str(node)+".txt | awk '{print $4}' >> launch_ids.txt", shell=True, executable='/bin/bash')
    elif args.run_type == 1:
        syscode = run()
        if syscode != 0:
            sys.exit(1)
    elif args.run_type == 2:
        unittest.main(argv=[''], verbosity=2, exit=False)
    else:
        assert False  # unrecognized run_type
