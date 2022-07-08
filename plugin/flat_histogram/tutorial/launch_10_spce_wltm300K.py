# This tutorial is similar to tutorial 5 spce, but for this low temperature simulation,
# we will split the simulation into two different nodes.
# The first node will have less particles but a higher number of sweeps required
# The second node will have dccb.

import sys
import subprocess
import argparse
import random
import unittest
import pathlib

# define parameters of a pure component SPCE simulation
# the remaining params depend on the node, and are thus defined in the run() function
params = {
    "cubic_box_length": 20, "fstprt": "/feasst/forcefield/spce.fstprt", "min_particles": 0,
    "temperature": 300, "max_particles": 296, "beta_mu": -15.24,
    "trials_per": 1e6, "hours_per_adjust": 0.01, "hours_per_checkpoint": 1, "seed": random.randrange(1e9), "num_hours": 5*24,
    "equilibration": 1e6, "num_nodes": 2, "procs_per_node": 32, "script": __file__, "dccb_cut": 0.9*3.165}
params["ewald_alpha"] = 5.6/params["cubic_box_length"]
R = 1.3806488E-23*6.02214129E+23 # J/mol/K
params["beta"] = 1./(params["temperature"]*R/1e3) # mol/kJ
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
def mc_spce(params, file_name):
    with open(file_name, "w") as myfile: myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_per_adjust} ln_prob_file spce_lnpin{node}.txt bounds_file spce_boundsn{node}.txt num_adjust_per_write 10
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha {window_alpha} min_size 2
Checkpoint file_name spce_checkpointn{node}.fst num_hours {hours_per_checkpoint} num_hours_terminate {num_hours_terminate}

# begin description of each MC clone
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt} physical_constants CODATA2010 \
    group0 oxygen oxygen_site_type 0
Potential VisitModel Ewald alpha {ewald_alpha} kmax_squared 38
Potential Model ModelTwoBodyFactory model0 LennardJones model1 ChargeScreened erfc_table_size 2e4 VisitModel VisitModelCutoffOuter
{ref_potential}
Potential Model ChargeScreenedIntra VisitModel VisitModelBond
Potential Model ChargeSelf
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 0.5 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 0.5 particle_type 0 tunable_param 0.5 tunable_target_acceptance 0.25
Log trials_per {trials_per} file_name spcen{node}s[sim_index].txt
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
{gce_trial}
RemoveAnalyze name Log
Log trials_per {trials_per} file_name spcen{node}s[sim_index].txt
Movie trials_per {trials_per} file_name spcen{node}s[sim_index].xyz
Tune trials_per_write {trials_per} file_name spce_tunen{node}s[sim_index].txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per} file_name spce_enn{node}s[sim_index].txt multistate true start_after_iteration 1
CriteriaUpdater trials_per 1e5
CriteriaWriter trials_per {trials_per} file_name spce_critn{node}s[sim_index].txt
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
class TestFlatHistogramSPCE(unittest.TestCase):
    def test(self):
        # compare the lnpi with the srsw
        import numpy as np
        import pandas as pd
        from pyfeasst import macrostate_distribution
        from pyfeasst import multistate_accumulator
        from pyfeasst import physical_constants
        num_block = 32
        # convert density of particles/A^3 to kg/m^3 of water
        rho_conv = 1e30/physical_constants.AvogadroConstant().value()*18.01528/1e3
        # convert pressure of kJ/mol/A^3 to bar
        pressure_conv = 1e30/physical_constants.AvogadroConstant().value()*1e3/1e5
        beta_mu_eq = list()
        rho_vapor = list()
        rho_liquid = list()
        en_vapor = list()
        en_liquid = list()
        pressure = list()
        for block in range(-1, num_block):
            print('block', block)
            if block == -1:
                ln_prob_header = 'ln_prob'
                energy_header = 'e_average'
            else:
                ln_prob_header = 'ln_prob' + str(block)
                energy_header = 'e_block' + str(block)
            lnpi = macrostate_distribution.splice_files(prefix='spce_lnpin', suffix='.txt',
                                                        ln_prob_header=ln_prob_header)
            lnpi.set_minimum_smoothing(60)
            energy = multistate_accumulator.splice_by_node(prefix='spce_enn', suffix='.txt', num_nodes=params['num_nodes'])
            lnpi.concat_dataframe(dataframe=energy, add_prefix='e_')
            delta_beta_mu = lnpi.equilibrium()
            beta_mu_eq.append(params['beta']*params['mu'] + delta_beta_mu)
            for index, phase in enumerate(lnpi.split()):
                n_gce = phase.average_macrostate()
                rho = rho_conv*n_gce/params['cubic_box_length']**3
                en = phase.ensemble_average(energy_header)/n_gce
                if index == 0:
                    pressure.append(-pressure_conv*phase.ln_prob()[0]/params['beta']/params['cubic_box_length']**3)
                    rho_vapor.append(rho)
                    en_vapor.append(en)
                else:
                    rho_liquid.append(rho)
                    en_liquid.append(en)
        data = pd.DataFrame(data={'rho_vapor': rho_vapor, 'rho_liquid': rho_liquid, 'pressure': pressure, 'en_vapor': en_vapor, 'en_liquid': en_liquid, 'beta_mu_eq': beta_mu_eq})
        data.to_csv('launch_10_spce_wltm300K.csv')
        for col in data.columns:
            print(col, data[col][0], '+/-', data[col][1:].std()/np.sqrt(len(data[col][1:])), data[col][1:].mean())

        # check equilibrium properties from https://www.nist.gov/mml/csd/chemical-informatics-group/sat-tmmc-liquid-vapor-coexistence-properties-long-range
        z_factor = 6
        self.assertLess(np.abs(data['rho_vapor'][0] - 7.373E-03), z_factor*np.sqrt((data['rho_vapor'][1:].std()/np.sqrt(num_block))**2 + 3.253E-05**2))
        self.assertLess(np.abs(data['rho_liquid'][0] - 9.981E+02), z_factor*np.sqrt((data['rho_liquid'][1:].std()/np.sqrt(num_block))**2 + 2.928E+00**2))
        self.assertLess(np.abs(data['pressure'][0] - 1.017E-02), z_factor*np.sqrt((data['pressure'][1:].std()/np.sqrt(num_block))**2 + 4.516E-05**2))
        self.assertLess(np.abs(data['beta_mu_eq'][0] - -1.524E+01), z_factor*np.sqrt((data['beta_mu_eq'][1:].std()/np.sqrt(num_block))**2 + 5.724E-02**2))

        # check lnpi
        lnpi = macrostate_distribution.splice_files(prefix='spce_lnpin', suffix='.txt')
        srsw = pd.read_csv('../test/data/colMatb0.400908lnz-15.24', skiprows=18, header=None, delim_whitespace=True)
        df = pd.concat([lnpi.dataframe(), srsw], axis=1)
        df['deltalnPI']=df[1]-df[1].shift(1)
        df.to_csv('spce_lnpi.csv')
        diverged=df[df.deltalnPI-df.delta_ln_prob > z_factor*df.delta_ln_prob_stdev]
        print(diverged)
        self.assertTrue(len(diverged) == 0)

# run the simulation and, if complete, analyze.
def run():
    if args.task == 0:
        file_name = "spce_launch"+str(params["node"])+".txt"
        splice_particles=180
        if params['node'] == 0:
            params["max_particles"]=splice_particles
            params["gce_trial"]="TrialTransfer weight 2 particle_type 0"
            params["ref_potential"]=""
            params["min_sweeps"]=5000
            params["window_alpha"]=1.1
        if params['node'] == 1:
            params["min_particles"]=splice_particles
            params["gce_trial"]="TrialGrowFile file_name spce_grow.txt"
            params["ref_potential"]="""RefPotential Model HardSphere group oxygen cutoff {dccb_cut} VisitModel VisitModelCell min_length {dccb_cut} cell_group oxygen""".format(**params)
            params["min_sweeps"]=50
            params["window_alpha"]=2.5
        mc_spce(params=params, file_name=file_name)
        syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > spce_launch"+str(params['node'])+".log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../../../build/bin/rst spce_checkpointn"+str(params['node'])+".fst", shell=True, executable='/bin/bash')
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
    else:
        assert False  # unrecognized run_type
