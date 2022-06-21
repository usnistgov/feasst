import sys
import subprocess
import argparse
import random
import unittest

# define parameters of a pure component NVT MC Lennard-Jones simulation
params = {
    "cubic_box_length": 8, "fstprt": "/feasst/forcefield/lj.fstprt", "beta": 1/1.5,
    "max_particles": 370, "min_particles": 0, "min_sweeps": 1e4, "mu": -2.352321,
    "trials_per": 1e6, "hours_per_adjust": 0.01, "hours_per_checkpoint": 1, "seed": random.randrange(1e9), "num_hours": 5*24,
    "equilibration": 1e6, "num_nodes": 1, "procs_per_node": 32, "script": __file__, "dccb_cut": 2**(1./6.)}
params["num_minutes"] = round(params["num_hours"]*60)
params["hours_per_adjust"] = params["hours_per_adjust"]*params["procs_per_node"]
params["hours_per_checkpoint"] = params["hours_per_checkpoint"]*params["procs_per_node"]
params["num_hours_terminate"] = 0.95*params["num_hours"]*params["procs_per_node"]
params["mu_init"] = params["mu"] + 1

# write fst script to run a single simulation
def mc_lj(params=params, file_name="launch.txt"):
    with open(file_name, "w") as myfile: myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_per_adjust} ln_prob_file lj_lnpi.txt bounds_file lj_bounds.txt num_adjust_per_write 10
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 2.5 min_size 2
Checkpoint file_name lj_checkpoint.fst num_hours {hours_per_checkpoint} num_hours_terminate {num_hours_terminate}

# begin description of each MC clone
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt}
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ConvertToRefPotential cutoff {dccb_cut} use_cell true
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
Log trials_per {trials_per} file_name lj[sim_index].txt
Tune
CheckEnergy trials_per {trials_per} tolerance 1e-8
#Checkpoint file_name lj_checkpoint[sim_index].fst num_hours {hours_per_checkpoint}

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Run until_num_particles [soft_macro_min]
RemoveTrial name TrialAdd
Run num_trials {equilibration}
RemoveModify name Tune

# gcmc tm production
ThermoParams beta {beta} chemical_potential {mu}
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} new_sweep 1 min_flatness 25 collect_flatness 20 min_collect_sweeps 20
TrialTransfer weight 2 particle_type 0 reference_index 0 num_steps 4
Tune trials_per_write {trials_per} file_name lj_tune[sim_index].txt multistate true
Movie trials_per {trials_per} file_name lj[sim_index].xyz
Energy trials_per_write {trials_per} file_name lj_en[sim_index].txt multistate true
CriteriaUpdater trials_per {trials_per}
CriteriaWriter trials_per {trials_per} file_name lj_crit[sim_index].txt
#Run until_criteria_complete true
""".format(**params))

# write slurm script to fill nodes with simulations
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

# after the simulation is complete, perform some tests
class TestFlatHistogramLJ(unittest.TestCase):
    def test(self):
        # analyze grand canonical ensemble average number of particles
        import numpy as np
        import pandas as pd
        lnpi=pd.read_csv('lj_lnpi.txt')
        self.assertAlmostEqual(310.4179421879679, (np.exp(lnpi["ln_prob"]) * lnpi["state"]).sum(), delta=0.25)

# run the simulation and, if complete, analyze.
def run():
    if args.task == 0:
        file_name = "lj_launch.txt"
        mc_lj(params, file_name=file_name)
        syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > lj_launch.log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../../../build/bin/rst lj_checkpoint.fst", shell=True, executable='/bin/bash')
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
