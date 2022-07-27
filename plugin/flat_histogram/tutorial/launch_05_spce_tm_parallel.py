import sys
import subprocess
import argparse
import random
import unittest
from pyfeasst import physical_constants

# define parameters of a pure component NVT MC SPCE simulation
params = {
    "cubic_box_length": 20, "fstprt": "/feasst/forcefield/spce.fstprt", "min_particles": 0,
    "temperature": 525, "max_particles": 265,  "min_sweeps": 400, "beta_mu": -8.14,
    #"temperature": 300, "max_particles": 296,  "min_sweeps": 200, "beta_mu": -15.24,
    "trials_per": 1e6, "hours_per_adjust": 0.01, "hours_per_checkpoint": 1, "seed": random.randrange(1e9), "num_hours": 5*24,
    "equilibration": 1e6, "num_nodes": 1, "procs_per_node": 32, "script": __file__, "dccb_cut": 0.9*3.165}
params["alpha"] = 5.6/params["cubic_box_length"]
params["beta"] = 1./(params["temperature"]*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
params["mu"] = params["beta_mu"]/params["beta"]
params["num_minutes"] = round(params["num_hours"]*60)
params["hours_per_adjust"] = params["hours_per_adjust"]*params["procs_per_node"]
params["hours_per_checkpoint"] = params["hours_per_checkpoint"]*params["procs_per_node"]
params["num_hours_terminate"] = 0.95*params["num_hours"]*params["procs_per_node"]
params["mu_init"] = -7
params["dccb_cut"] = params["cubic_box_length"]/int(params["cubic_box_length"]/params["dccb_cut"]) # maximize inside box

# write TrialGrowFile for SPCE
with open('spce_grow.txt', 'w') as f:
    f.write("""TrialGrowFile

particle_type 0 weight 2 transfer true site 0 num_steps 10 reference_index 0
bond true mobile_site 1 anchor_site 0 reference_index 0
angle true mobile_site 2 anchor_site 0 anchor_site2 1 reference_index 0
""")

# write fst script
def mc_spce(params=params, file_name="launch.txt"):
    with open(file_name, "w") as myfile: myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_per_adjust} ln_prob_file spce_lnpi.txt bounds_file spce_bounds.txt num_adjust_per_write 10
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 1.75 min_size 2
Checkpoint file_name spce_checkpoint.fst num_hours {hours_per_checkpoint} num_hours_terminate {num_hours_terminate}

# begin description of each MC clone
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt} physical_constants CODATA2010 \
    group0 oxygen oxygen_site_type 0
Potential VisitModel Ewald alpha {alpha} kmax_squared 38
Potential Model ModelTwoBodyFactory model0 LennardJones model1 ChargeScreened erfc_table_size 2e4 VisitModel VisitModelCutoffOuter
RefPotential Model HardSphere group oxygen cutoff {dccb_cut} VisitModel VisitModelCell min_length {dccb_cut} cell_group oxygen
Potential Model ChargeScreenedIntra VisitModel VisitModelBond
Potential Model ChargeSelf
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 0.5 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 0.5 particle_type 0 tunable_param 0.5 tunable_target_acceptance 0.25
Log trials_per {trials_per} file_name spce[sim_index].txt
Tune
CheckEnergy trials_per {trials_per} tolerance 1e-4

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Run until_num_particles [soft_macro_min]
RemoveTrial name TrialAdd
ThermoParams beta {beta} chemical_potential {mu}
Metropolis
Run num_trials {equilibration}
RemoveModify name Tune

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} new_sweep 1 min_flatness 25 collect_flatness 20 min_collect_sweeps 20
TrialGrowFile file_name spce_grow.txt
RemoveAnalyze name Log
Log trials_per {trials_per} file_name spce[sim_index].txt
Movie trials_per {trials_per} file_name spce[sim_index].xyz
Tune trials_per_write {trials_per} file_name spce_tune[sim_index].txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per} file_name spce_en[sim_index].txt multistate true start_after_iteration 1
CriteriaUpdater trials_per 1e5
CriteriaWriter trials_per {trials_per} file_name spce_crit[sim_index].txt
""".format(**params))

# write slurm script to fill node(s) with simulations
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
class TestFlatHistogramSPCE(unittest.TestCase):
    def test(self):
        # compare the lnpi with the srsw
        import numpy as np
        import pandas as pd
        df=pd.read_csv('spce_lnpi.txt')
        df=pd.concat([df, pd.read_csv('../test/data/stat_spce_525.csv')], axis=1)
        df['deltalnPI'] = df.lnPI - df.lnPI.shift(1)
        diverged=df[df.deltalnPI-df.delta_ln_prob > 6*df.delta_ln_prob_stdev]
        print(diverged)
        self.assertTrue(len(diverged) == 0)

# run the simulation and, if complete, analyze.
def run():
    if args.task == 0:
        file_name = "spce_launch.txt"
        mc_spce(params, file_name=file_name)
        syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > spce_launch.log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../../../build/bin/rst spce_checkpoint.fst", shell=True, executable='/bin/bash')
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
