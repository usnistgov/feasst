# This tutorial is similar to tutorial 10, but with a TraPPE alkane.

import sys
import subprocess
import argparse
import json
import random
import unittest
import pathlib
from pyfeasst import physical_constants

# define parameters of a pure component LJ simulation
#data=~/feasst/forcefield/n-butane.fstprt; temperature=350; max_particles=575; box_length=45; beta_mu=-6; cutoff=12; # https://www.nist.gov/mml/csd/chemical-informatics-research-group/sat-tmmc-liquid-vapor-coexistence-properties-trappe-ua-n
params = {
    #"cubic_box_length": 30, "fstprt": "/feasst/forcefield/ethane.fstprt", "temp_in_K": 300, 'max_particles': 225, 'betamu': -7.1126, 'num_sites': 2, 'cutoff': 14,
    #"cubic_box_length": 28, "fstprt": "/feasst/forcefield/propane.fstprt", "temp_in_K": 344, 'max_particles': 180, 'betamu': -7, 'num_sites': 3, 'cutoff': 14,

    # https://www.nist.gov/mml/csd/chemical-informatics-research-group/sat-tmmc-liquid-vapor-coexistence-properties-trappe-ua-n
    "cubic_box_length": 45, "fstprt": "/feasst/forcefield/n-butane.fstprt", "temp_in_K": 350, 'max_particles': 128, 'betamu': -6, 'num_sites': 4, 'cutoff': 12, 'molec_weight': 58.12,

    # https://mmlapps.nist.gov/srs/OCTANE/octane_sat.htm
    #"cubic_box_length": 35, "fstprt": "/feasst/forcefield/n-octane.fstprt", "temp_in_K": 560, 'max_particles': 128, 'betamu': -6, 'num_sites': 8, 'cutoff': 15,
    #"cubic_box_length": 36, "fstprt": "/feasst/forcefield/n-decane.fstprt", "temp_in_K": 560, 'max_particles': 12, 'betamu': -7, 'num_sites': 10, 'cutoff': 14,
    "dccb_cut": 4,
    "trials_per": 1e4, "hours_per_adjust": 10, "hours_per_checkpoint": 1, "seed": random.randrange(int(1e9)), "num_hours": 5*24,
    "equilibration": 1e5, "num_nodes": 1, "procs_per_node": 32, "script": __file__}
params["beta"] = 1./(params["temp_in_K"]*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
params['mu'] = params['betamu']/params['beta']
params["num_minutes"] = round(params["num_hours"]*60)
params["hours_per_adjust"] = params["hours_per_adjust"]*params["procs_per_node"]
params["hours_per_checkpoint"] = params["hours_per_checkpoint"]*params["procs_per_node"]
params["num_hours_terminate"] = 0.95*params["num_hours"]*params["procs_per_node"]
params["dccb_cut"] = params["cubic_box_length"]/int(params["cubic_box_length"]/params["dccb_cut"]) # maximize inside box
params["mu_init"]=10
params['last_site'] = params['num_sites'] - 1
def per_node_params():
    if params['num_nodes'] == 1:
        splice_particles=params['max_particles']
    elif params['num_nodes'] == 2:
        splice_particles=int(params['max_particles']*0.75)
    else:
        print('unrecognized processors per node:', params['procs_per_node'])
    if params['node'] == 0:
        params["min_particles"]=0
        params["max_particles"]=splice_particles
        params["min_sweeps"]=200
        params["window_alpha"]=2
        params["min_window_size"]=3
    elif params['node'] == 1:
        params["min_particles"]=splice_particles
        #params["max_particles"]=25
        params["min_sweeps"]=200
        params["window_alpha"]=1
        params["min_window_size"]=3
    else:
        assert False # unrecognized number of nodes

def write_partial(f, bond, angle, dihedral):
    if params['num_sites'] == 2:
        f.write(bond)
    elif params['num_sites'] == 3:
        f.write(angle)
    elif params['num_sites'] > 3:
        f.write(dihedral)
    else:
        print('unrecognized num_sites', params['num_sites'])
        assert False

# write TrialGrowFile to include grand canonical ensemble growth and canonica ensemble reptations
def write_grow_file(filename, gce):
    with open(filename, 'w') as f:
        f.write("TrialGrowFile\n\n")
        for inv in [True, False]:
            for trial_type in [0, 1, 2]: # 0: reptate, 1: full regrow, 2: partial regrow
                for site in range(params['num_sites']):
                    for i in range(4):
                        sign = -1
                        if (trial_type == 0 or trial_type == 2) and site != params['num_sites'] - 1:
                            sign = 1
                        params['site'+str(i)] = site + sign*i
                        if inv:
                            params['site'+str(i)] = params['num_sites'] - site - 1 - sign*i
                    bond = """bond true mobile_site {site0} anchor_site {site1} reference_index 0 num_steps 4 reference_index 0\n""".format(**params)
                    angle = """angle true mobile_site {site0} anchor_site {site1} anchor_site2 {site2} num_steps 4 reference_index 0\n""".format(**params)
                    dihedral = """dihedral true mobile_site {site0} anchor_site {site1} anchor_site2 {site2} anchor_site3 {site3} num_steps 4 reference_index 0\n""".format(**params)

                    # full regrowth insertion/deletion
                    if trial_type == 1 and gce:
                        if site == 0:
                            f.write("""particle_type 0 weight 2 transfer true site {site0} num_steps 4 reference_index 0\n""".format(**params))
                        elif site == 1:
                            f.write(bond)
                        elif site == 2:
                            f.write(angle)
                        else:
                            f.write(dihedral)

                    # reptation
                    elif trial_type == 0 and not gce:
                        if site == params['num_sites'] - 1:
                            write_partial(f, bond, angle, dihedral)
                        else:
                            if site == 0:
                                f.write("""particle_type 0 weight 2 """)
                            f.write("""reptate true mobile_site {site0} anchor_site {site1} num_steps 1 reference_index 0\n""".format(**params))

                    # partial regrow of the last site
                    if not gce and trial_type == 2:
                        if site == 0:
                            f.write("""particle_type 0 weight 2 """)
                            write_partial(f, bond, angle, dihedral)

                f.write("\n")

write_grow_file(filename="trappe_grow_canonical.txt", gce=False)
write_grow_file(filename="trappe_grow_grand_canonical.txt", gce=True)

# write fst script
def mc_trappe(params, file_name):
    with open(file_name, "w") as myfile: myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_per_adjust} ln_prob_file trappe_lnpin{node}.txt bounds_file trappe_boundsn{node}.txt num_adjust_per_write 10 min_window_size {min_window_size}
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha {window_alpha} min_size {min_window_size}
Checkpoint file_name trappe_checkpointn{node}.fst num_hours {hours_per_checkpoint} num_hours_terminate {num_hours_terminate}

# begin description of each MC clone
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt} cutoff {cutoff}
Potential Model LennardJones
Potential Model LennardJones VisitModel VisitModelIntra intra_cut 4
Potential VisitModel LongRangeCorrections
RefPotential Model LennardJones VisitModel VisitModelCell min_length {dccb_cut}
RefPotential Model LennardJones VisitModel VisitModelIntra intra_cut 4
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 0.5 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 0.25 particle_type 0 tunable_param 0.2 tunable_target_acceptance 0.25 pivot_site 0
TrialParticlePivot weight 0.25 particle_type 0 tunable_param 0.2 tunable_target_acceptance 0.25 pivot_site {last_site}
TrialGrowFile file_name trappe_grow_canonical.txt
Log trials_per {trials_per} file_name trappen{node}s[sim_index].txt
Tune
CheckEnergy trials_per {trials_per} tolerance 1e-4

# gcmc initialization and nvt equilibration
TrialGrowFile file_name trappe_grow_grand_canonical.txt
Run until_num_particles [soft_macro_min]
# Remove the 4 grand canonical trials for canonical equilibration, then add them back for production.
# Each transfer is two trials, two for each of forward and reverse order.
RemoveTrial index 7
RemoveTrial index 7
RemoveTrial index 7
RemoveTrial index 7
ThermoParams beta {beta} chemical_potential {mu}
Metropolis
Run num_trials {equilibration}
RemoveModify name Tune

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} new_sweep 1 min_flatness 25 collect_flatness 20 min_collect_sweeps 20
TrialGrowFile file_name trappe_grow_grand_canonical.txt
RemoveAnalyze name Log
Log trials_per {trials_per} file_name trappen{node}s[sim_index].txt
#Movie trials_per {trials_per} file_name trappen{node}s[sim_index].xyz
Tune trials_per_write {trials_per} file_name trappe_tunen{node}s[sim_index].txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per} file_name trappe_enn{node}s[sim_index].txt multistate true start_after_iteration 1
CriteriaUpdater trials_per 1e5
CriteriaWriter trials_per {trials_per} file_name trappe_critn{node}s[sim_index].txt
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
if args.params != "":
    with open(args.params) as jsonfile:
        params = json.loads(json.load(jsonfile))

# after the simulation is complete, perform some analysis
class TestFlatHistogramLJ(unittest.TestCase):
    def test(self):
        from pyfeasst import macrostate_distribution
        from pyfeasst import physical_constants
        lnp = macrostate_distribution.splice_collection_matrix(prefix='trappe_critn0s', suffix='.txt', use_soft=True)
        lnp.equilibrium()
        vapor, liquid = lnp.split()
        volume = params['cubic_box_length']**3
        na = physical_constants.AvogadroConstant().value()
        dens_conv = 1./volume/na*params['molec_weight']/1e3*1e30 # convert from N/V units of molecules/A^3 to kg/m
        # https://www.nist.gov/mml/csd/chemical-informatics-research-group/sat-tmmc-liquid-vapor-coexistence-properties-trappe-ua-n
        self.assertAlmostEqual(30.6, vapor.average_macrostate()*dens_conv, delta=1)
        self.assertAlmostEqual(508, liquid.average_macrostate()*dens_conv, delta=30)
        #lnp.plot(show=True)

# run the simulation and, if complete, analyze.
def run():
    if args.task == 0:
        file_name = "trappe_launch"+str(params["node"])+".txt"
        mc_trappe(params=params, file_name=file_name)
        syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > trappe_launch"+str(params['node'])+".log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../../../build/bin/rst trappe_checkpointn"+str(params['node'])+".fst", shell=True, executable='/bin/bash')
    if syscode == 0 and params['node'] == 1:
        unittest.main(argv=[''], verbosity=2, exit=False)
    return syscode

if __name__ == "__main__":
    for node in range(params['num_nodes']):
        print('node', node)
        params["node"] = node
        params['params'] = 'trappe_params' + str(node) + '.txt'
        per_node_params()
    if args.run_type == 0:
            with open(params['params'], 'w') as jsonfile:
                json.dump(json.dumps(params), jsonfile)
            slurm_file = 'trappe_slurm'+str(node)+'.txt'
            slurm_queue(slurm_file)
            subprocess.call("sbatch --array=0-1%1 " + slurm_file + " | awk '{print $4}' >> launch_ids.txt", shell=True, executable='/bin/bash')
    elif args.run_type == 1:
        syscode = run()
        if syscode != 0:
            sys.exit(1)
    elif args.run_type == 2:
        unittest.main(argv=[''], verbosity=2, exit=False)
    else:
        assert False  # unrecognized run_type
