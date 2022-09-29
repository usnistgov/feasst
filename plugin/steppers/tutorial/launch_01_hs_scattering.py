import sys
import subprocess
import argparse
import random
import unittest

# define parameters of a pure component NVT MC hard sphere simulation
params = {
    "cubic_box_length": 8, "fstprt": "/feasst/forcefield/atom.fstprt",
    "num_particles": 128,
    "trials_per": 1e6, "hours_per_checkpoint": 1, "seed": random.randrange(int(1e9)), "num_hours": 5*24,
    "equilibration": 1e5, "production": 1e6, "num_nodes": 1, "procs_per_node": 1, "script": __file__}
params["num_minutes"] = round(params["num_hours"]*60)
params["num_hours_terminate"] = 0.95*params["num_hours"]*params["procs_per_node"]

# write fst script to run a single simulation
def mc_hs(params=params, file_name="launch.txt"):
    with open(file_name, "w") as myfile: myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt}
Potential Model HardSphere VisitModel VisitModelCell min_length 1
ThermoParams beta 1 chemical_potential 1
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
Log trials_per {trials_per} file_name hs.csv
Tune
CheckEnergy trials_per {trials_per} tolerance 1e-8
Checkpoint file_name hs.fst num_hours_terminate {num_hours_terminate}

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Run until_num_particles {num_particles}
RemoveTrial name TrialAdd
Run num_trials {equilibration}

# nvt production
Movie trials_per {trials_per} file_name hs.xyz
PairDistribution trials_per_update 1000 trials_per_write {trials_per} dr 0.025 file_name hs_gr.csv
Scattering trials_per_update 100 trials_per_write {trials_per} num_frequency 10 file_name hs_iq.csv
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
        from pyfeasst import scattering
        gr=pd.read_csv('hs_gr.csv', comment="#")
        iq=pd.read_csv('hs_iq.csv', comment="#")
        grp = iq.groupby('q', as_index=False)
        self.assertAlmostEqual(gr['g0-0'][45], 1.2829, delta=0.05)
        self.assertAlmostEqual(iq['i'][3810], 5.72894, delta=0.4)
        self.assertAlmostEqual(iq['i'][0]/iq['p0'][0]**2, 1, delta=0.075)

        # scale the gr closer to one at the tail by dividing by the average of the last 5%
        # A fourier transform of this scaled gr will result in a smoother low-q
        gr_scaled = gr['g0-0']/np.average(gr['g0-0'][int(0.95*len(gr['r'])):])
        import matplotlib.pyplot as plt
        plt.plot(gr['r'], gr['g0-0'], label='gr')
        plt.scatter(iq['q'], iq['i']/iq['p0']**2, label='sq_all')
        plt.plot(grp.mean()['q'], grp.mean()['i']/grp.mean()['p0']**2, label='sq_av')
        qs = np.arange(2*np.pi/8, 10, 0.01)
        sq=list(); sqs=list()
        number_density = params['num_particles']/params['cubic_box_length']**3
        for q in qs:
            sqs.append(scattering.structure_factor(gr['r'], gr_scaled, frequency=q, number_density=number_density))
            sq.append(scattering.structure_factor(gr['r'], gr['g0-0'], frequency=q, number_density=number_density))
        plt.plot(qs, sq, label='sq_from_gr')
        plt.plot(qs, sqs, label='sq_from_scaled_gr')
        plt.legend()
        #plt.show()

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
    else:
        assert False  # unrecognized run_type
