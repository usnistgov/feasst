# This tutorial is similar to tutorial 4 lj, but for this low temperature simulation,
# we will split the simulation into two different nodes.
# The first node will have less particles but a higher number of sweeps required
# The second node will have dccb but not avb.

import sys
import subprocess
import argparse
import json
import random
import unittest
import pathlib

# define parameters of a pure component LJ simulation
params = {
    "cubic_box_length": 8, "fstprt": "/feasst/forcefield/lj.fstprt", "dccb_cut": 1, "beta": 1/0.7, "mu": -4.1603632,
    "trials_per": 1e6, "hours_per_adjust": 0.01, "hours_per_checkpoint": 1, "seed": random.randrange(int(1e9)), "num_hours": 5*24,
    "equilibration": 1e6, "num_nodes": 2, "procs_per_node": 32, "script": __file__}
params["num_minutes"] = round(params["num_hours"]*60)
params["hours_per_adjust"] = params["hours_per_adjust"]*params["procs_per_node"]
params["hours_per_checkpoint"] = params["hours_per_checkpoint"]*params["procs_per_node"]
params["num_hours_terminate"] = 0.95*params["num_hours"]*params["procs_per_node"]
params["dccb_cut"] = params["cubic_box_length"]/int(params["cubic_box_length"]/params["dccb_cut"]) # maximize inside box
params["mu_init"]=10
def per_node_params():
    splice_particles=375
    if params['node'] == 0:
        params["min_particles"]=0
        params["max_particles"]=splice_particles
        params["gce_trial"]="TrialTransfer weight 2 particle_type 0\nTrialTransferAVB weight 0.2 particle_type 0"
        params["lj_potential"]="Potential EnergyMap EnergyMapNeighborCriteria neighbor_index 0 Model LennardJones"
        params["ref_potential"]=""
        params["avb_trials"]="TrialAVB2 weight 0.1 particle_type 0\nTrialAVB4 weight 0.1 particle_type 0"
        params["min_sweeps"]=2000
        params["window_alpha"]=2
        params["min_window_size"]=5
    elif params['node'] == 1:
        params["min_particles"]=splice_particles
        params["max_particles"]=475
        params["gce_trial"] = "TrialTransfer weight 2 particle_type 0 reference_index 0 num_steps 10"
        params["lj_potential"]="Potential Model LennardJones"
        params["ref_potential"]="""RefPotential Model LennardJones cutoff {dccb_cut} VisitModel VisitModelCell min_length {dccb_cut}""".format(**params)
        params["avb_trials"]=""
        params["min_sweeps"]=200
        params["window_alpha"]=1
        params["min_window_size"]=3
    else:
        assert False # unrecognized number of nodes

# write fst script
def mc_lj(params, file_name):
    with open(file_name, "w") as myfile: myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_per_adjust} ln_prob_file lj_lnpin{node}.txt bounds_file lj_boundsn{node}.txt num_adjust_per_write 10 min_window_size {min_window_size}
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha {window_alpha} min_size {min_window_size}
Checkpoint file_name lj_checkpointn{node}.fst num_hours {hours_per_checkpoint} num_hours_terminate {num_hours_terminate}

# begin description of each MC clone
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt}
NeighborCriteria maximum_distance 1.375 minimum_distance 0.9
{lj_potential}
{ref_potential}
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
{avb_trials}
Log trials_per_write {trials_per} file_name ljn{node}s[sim_index].txt
Tune
CheckEnergy trials_per_update {trials_per} tolerance 1e-4

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
{gce_trial}
RemoveAnalyze name Log
Log trials_per_write {trials_per} file_name ljn{node}s[sim_index].txt
Movie trials_per_write {trials_per} file_name ljn{node}s[sim_index].xyz
Tune trials_per_write {trials_per} file_name lj_tunen{node}s[sim_index].txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per} file_name lj_enn{node}s[sim_index].txt multistate true start_after_iteration 1
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {trials_per} file_name lj_critn{node}s[sim_index].txt
""".format(**params))

# write slurm script
def slurm_queue(file_name):
    with open(file_name, "w") as myfile: myfile.write("""#!/bin/bash
#SBATCH -n {procs_per_node} -N 1 -t {num_minutes}:00 -o hostname_%j.out -e hostname_%j.out
echo "Running {script} ID $SLURM_JOB_ID on $(hostname) at $(date) in $PWD"
cd $PWD
export OMP_NUM_THREADS={procs_per_node}
python {script} --run_type 1 --task $SLURM_ARRAY_TASK_ID --params {params}
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
parser.add_argument('--params', type=str, default="", help="file name of the params file.")
args = parser.parse_args()
if args.params is not "":
    with open(args.params) as jsonfile:
        params = json.loads(json.load(jsonfile))

# after the simulation is complete, perform some analysis
class TestFlatHistogramLJ(unittest.TestCase):
    def test(self):
        # compare the lnpi with the srsw
        import numpy as np
        import pandas as pd
        from pyfeasst import macrostate_distribution
        from pyfeasst import multistate_accumulator
        num_block = 32
        beta_mu_eq = list()
        rho_vapor = list()
        rho_liquid = list()
        en_vapor = list()
        en_liquid = list()
        pressure = list()
        for block in range(-1, num_block):
            if block == -1:
                ln_prob_header = 'ln_prob'
                energy_header = 'e_average'
            else:
                ln_prob_header = 'ln_prob' + str(block)
                energy_header = 'e_block' + str(block)
            lnpi = macrostate_distribution.splice_files(prefix='lj_lnpin', suffix='.txt',
                                                        ln_prob_header=ln_prob_header)
            lnpi.set_minimum_smoothing(30)
            energy = multistate_accumulator.splice_by_node(prefix='lj_enn', suffix='.txt', num_nodes=params['num_nodes'])
            lnpi.concat_dataframe(dataframe=energy, add_prefix='e_')
            delta_beta_mu = lnpi.equilibrium()
            beta_mu_eq.append(params['beta']*params['mu'] + delta_beta_mu)
            for index, phase in enumerate(lnpi.split()):
                n_gce = phase.average_macrostate()
                rho = n_gce/params['cubic_box_length']**3
                en = phase.ensemble_average(energy_header)/n_gce
                if index == 0:
                    pressure.append(-phase.ln_prob()[0]/params['beta']/params['cubic_box_length']**3)
                    rho_vapor.append(rho)
                    en_vapor.append(en)
                else:
                    rho_liquid.append(rho)
                    en_liquid.append(en)
        data = pd.DataFrame(data={'rho_vapor': rho_vapor, 'rho_liquid': rho_liquid, 'pressure': pressure, 'en_vapor': en_vapor, 'en_liquid': en_liquid, 'beta_mu_eq': beta_mu_eq})
        data.to_csv('launch_10_lj_wltm_t0.7.csv')
        for col in data.columns:
            print(col, data[col][0], '+/-', data[col][1:].std()/np.sqrt(len(data[col][1:])), data[col][1:].mean())

        # check equilibrium properties from https://www.nist.gov/mml/csd/chemical-informatics-group/sat-tmmc-liquid-vapor-coexistence-properties-long-range
        z_factor = 6
        self.assertLess(np.abs(data['rho_vapor'][0] - 1.996E-03), z_factor*np.sqrt((data['rho_vapor'][1:].std()/np.sqrt(num_block))**2 + 1.422E-05**2))
        self.assertLess(np.abs(data['rho_liquid'][0] - 8.437E-01), z_factor*np.sqrt((data['rho_liquid'][1:].std()/np.sqrt(num_block))**2 + 2.49E-04**2))
        self.assertLess(np.abs(data['pressure'][0] - 1.370E-03), z_factor*np.sqrt((data['pressure'][1:].std()/np.sqrt(num_block))**2 + 9.507E-07**2))
        self.assertLess(np.abs(data['en_vapor'][0] - -2.500E-02), z_factor*np.sqrt((data['en_vapor'][1:].std()/np.sqrt(num_block))**2 + 3.323E-05**2))
        self.assertLess(np.abs(data['en_liquid'][0] - -6.106E+00), z_factor*np.sqrt((data['en_liquid'][1:].std()/np.sqrt(num_block))**2 + 1.832E-03**2))
        self.assertLess(np.abs(data['beta_mu_eq'][0] - -6.257E+00), z_factor*np.sqrt((data['beta_mu_eq'][1:].std()/np.sqrt(num_block))**2 + 4.639E-04**2))

        # check lnpi
        lnpi = macrostate_distribution.splice_files(prefix='lj_lnpin', suffix='.txt')
        df=pd.concat([lnpi.dataframe(), pd.read_csv('../test/data/stat070.csv')], axis=1)
        df['deltalnPI']=df.lnPI-df.lnPI.shift(1)
        df.to_csv('lj_lnpi.csv')
        diverged=df[df.deltalnPI-df.delta_ln_prob > z_factor*df.delta_ln_prob_stdev]
        print(diverged)
        self.assertTrue(len(diverged) == 0)

# run the simulation and, if complete, analyze.
def run():
    if args.task == 0:
        file_name = "lj_launch"+str(params["node"])+".txt"
        mc_lj(params=params, file_name=file_name)
        syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > lj_launch"+str(params['node'])+".log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../../../build/bin/rst lj_checkpointn"+str(params['node'])+".fst", shell=True, executable='/bin/bash')
    if syscode == 0 and params['node'] == 1:
        unittest.main(argv=[''], verbosity=2, exit=False)
    return syscode

if __name__ == "__main__":
    if args.run_type == 0:
        for node in range(params['num_nodes']):
            params["node"] = node
            params['params'] = 'lj_params'+str(node)
            per_node_params()
            with open(params['params'], 'w') as jsonfile:
                json.dump(json.dumps(params), jsonfile)
            slurm_file = 'lj_slurm'+str(node)+'.txt'
            slurm_queue(slurm_file)
            subprocess.call("sbatch --array=0-1%1 " + slurm_file + " | awk '{print $4}' >> launch_ids.txt", shell=True, executable='/bin/bash')
    elif args.run_type == 1:
        syscode = run()
        if syscode != 0:
            sys.exit(1)
    else:
        assert False  # unrecognized run_type
