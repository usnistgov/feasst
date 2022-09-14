import sys
import subprocess
import argparse
import math
import random
import unittest

# define parameters of a pure component NVT MC RPM expanded GCE simulation.
# See https://doi.org/10.1063/1.5123683
params = {
    "cubic_box_length": 12, "beta": 1./0.047899460618081,
    "plus": "/feasst/plugin/charge/forcefield/rpm_plus.fstprt",
    "minus": "/feasst/plugin/charge/forcefield/rpm_minus.fstprt",
    "max_particles": 20, "min_particles": 0, "min_sweeps": 1e3, "beta_mu": -13.94,
    "trials_per": 1e6, "hours_per_adjust": 0.01, "hours_per_checkpoint": 1, "seed": random.randrange(int(1e9)), "num_hours": 5*24,
    "equilibration": 1e6, "num_nodes": 1, "procs_per_node": 4, "script": __file__, "dccb_cut": 2**(1./6.)}
params["alpha"] = 6.87098396396261/params["cubic_box_length"]
params["mu"] = params["beta_mu"]/params["beta"]
params["charge_plus"] = 1./math.sqrt(1.602176634E-19**2/(4*math.pi*8.8541878128E-12*1e3/1e10/6.02214076E+23))
params["charge_minus"] = -params["charge_plus"]
params["num_minutes"] = round(params["num_hours"]*60)
params["hours_per_adjust"] = params["hours_per_adjust"]*params["procs_per_node"]
params["hours_per_checkpoint"] = params["hours_per_checkpoint"]*params["procs_per_node"]
params["num_hours_terminate"] = 0.95*params["num_hours"]*params["procs_per_node"]

# write fst script to run a single simulation
def mc_rpm(params=params, file_name="launch.txt"):
    with open(file_name, "w") as myfile: myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice min_window_size 2 hours_per {hours_per_adjust} ln_prob_file lnpi_rpm.txt bounds_file bounds_rpm.txt num_adjust_per_write 10
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 2.5
Checkpoint file_name checkpoint.fst num_hours {hours_per_checkpoint} num_hours_terminate {num_hours_terminate}

# begin description of each MC clone
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {plus} particle_type1 {minus} cutoff 4.891304347826090 charge0 {charge_plus} charge1 {charge_minus}
Potential VisitModel Ewald alpha {alpha} kmax_squared 38
Potential Model ModelTwoBodyFactory model0 HardSphere model1 ChargeScreened table_size 1e6
Potential Model ChargeSelf
ConvertToRefPotential potential_index 1 cutoff {dccb_cut} use_cell true
ThermoParams beta {beta} chemical_potential0 {mu} chemical_potential1 {mu}
Metropolis Constraint AEqualB extra_A 1
TrialTranslate weight 0.5 tunable_param 0.2 tunable_target_acceptance 0.25
Log trials_per {trials_per} file_name rpm[sim_index].txt
Tune
CheckEnergy trials_per {trials_per} tolerance 1e-4
#Checkpoint file_name checkpoint[sim_index].fst num_hours {hours_per_checkpoint}

# gcmc initialization and nvt equilibration
TrialAdd weight 1 particle_type 0 reference_index 0 num_steps 8
TrialAdd weight 1 particle_type 1 reference_index 0 num_steps 8
Run until_num_particles [soft_macro_min]
RemoveTrial name TrialAdd
RemoveTrial name TrialAdd
Run num_trials {equilibration}
RemoveModify name Tune

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
  Bias TransitionMatrix min_sweeps {min_sweeps} new_sweep 1 \
  Constraint AEqualB extra_A 1
#Bias WLTM min_sweeps {min_sweeps} new_sweep 1 min_flatness 25 collect_flatness 20
TrialTransfer weight 1 particle_type 0 reference_index 0 num_steps 8
TrialTransfer weight 1 particle_type 1 reference_index 0 num_steps 8
RemoveAnalyze name Log
Log trials_per {trials_per} file_name rpm[sim_index].txt
Tune trials_per_write {trials_per} file_name tune_rpm[sim_index].txt multistate true
Movie trials_per {trials_per} file_name rpm[sim_index].xyz
Energy trials_per_write {trials_per} file_name en_rpm[sim_index].txt multistate true
CriteriaUpdater trials_per {trials_per}
CriteriaWriter trials_per {trials_per} file_name crit_rpm[sim_index].txt
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
class TestFlatHistogramRPM(unittest.TestCase):
    def test(self):
        # compare the lnpi and energy with some known values
        import numpy as np
        import pandas as pd
        lnpi=pd.read_csv('lnpi_rpm.txt')
        lnpi=lnpi[:5] # cut down to five rows
        lnpi=lnpi[::2] # drop odd rows for neutral-only macrostates
        lnpi['ln_prob'] -= np.log(sum(np.exp(lnpi['ln_prob'])))  # renormalize
        lnpi['ln_prob_prev'] = [-1.2994315780357, -1.08646312498868, -0.941850889679828]
        lnpi['ln_prob_prev_stdev'] = [0.07, 0.05, 0.05]
        diverged=lnpi[lnpi.ln_prob-lnpi.ln_prob_prev > 5*lnpi.ln_prob_prev_stdev]
        self.assertTrue(len(diverged) == 0)
        en=pd.read_csv('en_rpm0.txt')
        en=en[:5]
        en['prev'] = [0, -0.115474, -0.939408, -1.32485, -2.02625]
        en['prev_stdev'] = [1e-14, 1e-6, 0.02, 0.03, 0.04]
        diverged=en[en.average-en.prev > 10*en.prev_stdev]
        self.assertTrue(len(diverged) == 0)


# run the simulation and, if complete, analyze.
def run():
    if args.task == 0:
        file_name = "launch_rpm.txt"
        mc_rpm(params, file_name=file_name)
        syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > launch_rpm.log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../../../build/bin/rst checkpoint.fst", shell=True, executable='/bin/bash')
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
        assert False  # unrecognized run_type
