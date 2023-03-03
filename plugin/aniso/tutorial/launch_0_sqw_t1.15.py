"""
This tutorial is similar to flat histogram tutorial 10, but for a square well
with anisotropic "texture" from random noise.
This is more of a mechanical test of the anisotropic potentials reproducing
an isotropic potential:
https://www.nist.gov/mml/csd/chemical-informatics-group/sat-tmmc-liquid-vapor-coexistence-properties-square-well-fluid
"""

import sys
import subprocess
import argparse
import random
import unittest
import numpy as np

params = {
    "cubic_box_length": 9, "fstprt": "/feasst/plugin/aniso/forcefield/aniso_tabular.fstprt", "dccb_cut": 1, "beta": 1/1.15, "mu": -3,
    "max_particles": 475, "min_particles": 0, "min_sweeps": 1e4,
    "trials_per": 1e6, "hours_per_adjust": 0.01, "hours_per_checkpoint": 1, "seed": random.randrange(int(1e9)), "num_hours": 5*24,
    "equilibration": 1e6, "num_nodes": 1, "procs_per_node": 32, "script": __file__,
    "min_window_size": 5, 'cutoff': 1.5}
params["num_minutes"] = round(params["num_hours"]*60)
params["hours_per_adjust"] = params["hours_per_adjust"]*params["procs_per_node"]
params["hours_per_checkpoint"] = params["hours_per_checkpoint"]*params["procs_per_node"]
params["num_hours_terminate"] = 0.95*params["num_hours"]*params["procs_per_node"]
params["mu_init"] = params["mu"] + 1

# generate the tabular potential
def generate_table(num_orientations_per_pi, num_z, file_name, gamma=1., delta=0.5):
    dt = np.pi/num_orientations_per_pi
    assert num_z > 1
    dz = 1./(num_z - 1)
    num_orientations = 0
    num_elements = 0
    fl = open(file_name, 'w')
    fl.write('site_types 1 0\n' +
             'num_orientations_per_pi ' + str(num_orientations_per_pi) + '\n' +
             'gamma ' + str(gamma) + '\n' +
             'delta ' + str(delta) + '\n' +
             'num_z ' + str(num_z) + '\n' +
             'smoothing_distance -1\n')
    #for s1 in np.arange(0, 2*np.pi + dt/2, dt): #theta
    for s1 in np.arange(0, np.pi + dt/2, dt): # if i==j, avoid x < 0 for i/j swap symmetry
        for s2 in np.arange(0, np.pi + dt/2, dt): #phi
            for e1 in np.arange(-np.pi, np.pi + dt/2, dt):
                for e2 in np.arange(0, np.pi + dt/2, dt):
                    for e3 in np.arange(-np.pi, np.pi + dt/2, dt):
                        num_orientations += 1
                        # hard code square well for testing
                        rh = 1.
                        fl.write(str(rh) + " ")
                        for z in np.arange(0, 1. + dz/2, dz):
                            num_elements += 1
                            # hard code square well for testing
                            energy = -1. + 0.0001*random.uniform(-1., 1.)
                            fl.write(str(energy) + " ")
                        fl.write("\n")
    fl.close()
    print('num_orientations', num_orientations)
    print('num_elements', num_elements)

# write fst script to run a single simulation
def mc_sqw(params=params, file_name="launch.txt"):
    generate_table(num_orientations_per_pi=3, num_z=2, file_name='dat.txt')
    with open(file_name, "w") as myfile: myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_per_adjust} ln_prob_file sqw_lnpi.txt bounds_file sqw_bounds.txt num_adjust_per_write 10 min_window_size {min_window_size}
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 2.25 min_size {min_window_size}
Checkpoint file_name sqw_checkpoint.fst num_hours {hours_per_checkpoint} num_hours_terminate {num_hours_terminate}

# begin description of each MC clone
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt} cutoff {cutoff}
Potential Model TwoBodyTable VisitModelInner VisitModelInnerTable table_file dat.txt
RefPotential Model HardSphere cutoff {dccb_cut} VisitModel VisitModelCell min_length {dccb_cut}
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
TrialRotate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
Log trials_per {trials_per} file_name sqw[sim_index].txt
Tune
CheckEnergy trials_per {trials_per} tolerance 1e-8
#Checkpoint file_name sqw_checkpoint[sim_index].fst num_hours {hours_per_checkpoint}

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
Movie trials_per {trials_per} file_name sqw[sim_index].xyz
Tune trials_per_write {trials_per} file_name sqw_tune[sim_index].txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per} file_name sqw_en[sim_index].txt multistate true start_after_iteration 1
CriteriaUpdater trials_per {trials_per}
CriteriaWriter trials_per {trials_per} file_name sqw_crit[sim_index].txt
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
parser.add_argument('--run_type', '-r', type=int, default=0, help="0: submit batch to scheduler, 1: run batch on host, 2: test")
parser.add_argument('--task', type=int, default=0, help="input by slurm scheduler. If >0, restart from checkpoint.")
args = parser.parse_args()

# after the simulation is complete, perform some tests
class TestFlatHistogramSQWT(unittest.TestCase):
    def test(self):
        # analyze grand canonical ensemble average number of particles
        import numpy as np
        import pandas as pd
        from pyfeasst import macrostate_distribution
        lnpi = macrostate_distribution.MacrostateDistribution(file_name='sqw_lnpi.txt')
        rw = lnpi.equilibrium()
        self.assertAlmostEqual(params['beta']*params['mu'] + rw, -3.194, delta=1e-3)
        #lnpi.plot(show=True)
        vap, liq = lnpi.split()
        self.assertAlmostEqual(vap.average_macrostate()/params['cubic_box_length']**3, 9.723E-02, delta=1e-3)
        self.assertAlmostEqual(liq.average_macrostate()/params['cubic_box_length']**3, 5.384E-01, delta=1e-3)

# run the simulation and, if complete, analyze.
def run():
    if args.task == 0:
        file_name = "sqw_launch.txt"
        mc_sqw(params, file_name=file_name)
        syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > sqw_launch.log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../../../build/bin/rst sqw_checkpoint.fst", shell=True, executable='/bin/bash')
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
