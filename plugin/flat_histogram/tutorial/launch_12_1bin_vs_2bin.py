"""
Comparison of the efficiency of 1 vs 2 bin flat histogram simulations, as described in
https://doi.org/10.1021/acs.jpcb.3c00613
"""

import sys
import subprocess
import argparse
import random
import unittest
import pathlib
import numpy as np
from multiprocessing import Pool
from pyfeasst import cd

# define parameters of a pure component NVT MC Lennard-Jones simulation
params = {
    "cubic_side_length": 8, "fstprt": "/feasst/particle/lj.fstprt", "beta": 1/0.7, "mu": -4.1603632,
    #"cubic_side_length": 8, "fstprt": "/feasst/particle/lj.fstprt", "beta": 1/1.5, "mu": -2.352321,
    "min_particles": 100, "min_sweeps": 1e9,
    "trials_per": 1e6, "hours_per_adjust": 0.01, "hours_per_checkpoint": 1, "seed": random.randrange(1e9), "num_hours": 0.5,
    "equilibration": 1e6, "num_nodes": 1, "procs_per_sim": 1, "procs_per_node": 4, "script": __file__, "dccb_cut": 2**(1./6.),
    "dir": str(pathlib.Path(__file__).parent.resolve())}
params["max_particles"] = params["min_particles"] + 1
params["num_minutes"] = round(params["num_hours"]*60)
params["hours_per_adjust"] = params["hours_per_adjust"]*params["procs_per_sim"]
params["hours_per_checkpoint"] = params["hours_per_checkpoint"]*params["procs_per_sim"]
params["num_hours_terminate"] = 0.95*params["num_hours"]*params["procs_per_sim"]
params["mu_init"] = params["mu"] + 1

# write fst script to run a single simulation
def mc_lj(params=params, script_file="launch.txt"):
    with open(script_file, "w") as myfile: myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_per_adjust} ln_prob_file lj_lnpin{min_particles}s{sim}.txt ln_prob_file_append true bounds_file lj_boundsn{min_particles}s{sim}.txt num_adjust_per_write 1 min_window_size 1
{window_custom}
Checkpoint checkpoint_file lj_checkpointn{min_particles}s{sim}.fst num_hours {hours_per_checkpoint} num_hours_terminate {num_hours_terminate}

# begin description of each MC clone
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt}
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
#ConvertToRefPotential cutoff {dccb_cut} use_cell true
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
Log trials_per_write {trials_per} output_file ljn{min_particles}s{sim}_[sim_index].txt
Tune
CheckEnergy trials_per_update {trials_per} tolerance 1e-8

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Run until_num_particles [soft_macro_min]
RemoveTrial name TrialAdd
Run num_trials {equilibration}
RemoveModify name Tune

# gcmc tm production
ThermoParams beta {beta} chemical_potential {mu}
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias TransitionMatrix min_sweeps {min_sweeps} new_sweep 1
#Bias WLTM min_sweeps {min_sweeps} new_sweep 1 min_flatness 25 collect_flatness 20 min_collect_sweeps 20
TrialTransfer weight 2 particle_type 0
#TrialTransfer weight 2 particle_type 0 reference_index 0 num_steps 4
Movie trials_per_write {trials_per} output_file ljn{min_particles}s{sim}_[sim_index].xyz
Tune trials_per_write {trials_per} output_file lj_tunen{min_particles}s{sim}_[sim_index].txt multistate true
Energy trials_per_write {trials_per} output_file lj_enn{min_particles}s{sim}_[sim_index].txt multistate true append true
CPUTime trials_per_write {trials_per} output_file lj_cpun{min_particles}s{sim}_[sim_index].txt append true
CriteriaUpdater trials_per_update {trials_per}
CriteriaWriter trials_per_write {trials_per} output_file lj_critn{min_particles}s{sim}_[sim_index].txt
#Run until_criteria_complete true
""".format(**params))

# write slurm script to fill nodes with simulations
def slurm_queue(min_particles):
    params["min_particles"] = min_particles
    with open("slurm.txt", "w") as myfile: myfile.write("""#!/bin/bash
#SBATCH -n {procs_per_node} -N {num_nodes} -t {num_minutes}:00 -o hostname_%j.out -e hostname_%j.out
echo "Running {script} ID $SLURM_JOB_ID on $(hostname) at $(date) in $PWD"
original_dir=$PWD; echo $original_dir
wrk=/wrk/$LOGNAME/$SLURM_JOB_ID/; mkdir -p $wrk; cd $wrk; echo "wrk:$wrk"
rsync -au $original_dir/* .; rm hostname_*
ls
export OMP_NUM_THREADS={procs_per_node}
python {script} --run_type 0 --dir $original_dir --min_particles {min_particles}
#python {script} --run_type 0 --task $SLURM_ARRAY_TASK_ID --dir $original_dir --min_particles {min_particles}
echo "sync"
rsync -au . $original_dir/
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
parser.add_argument('--run_type', '-r', type=int, default=0, help="0: run, 1: submit to queue")
#parser.add_argument('--task', type=int, default=0, help="input by slurm scheduler. If >0, restart from checkpoint.")
parser.add_argument('--dir', type=str, default="", help="orginal working directory that the script was launched.")
parser.add_argument('--min_particles', type=int, default=100, help="minimum number of particles")
args = parser.parse_args()
params['min_particles'] = args.min_particles
params["max_particles"] = params["min_particles"] + 1
if args.dir:
    params['dir']=args.dir

def dists_to_time_std(dists):
    time = list()
    std = list()
    for dist in dists:
        time.append(dist[1]['cpu_hours'])
        std.append(dist[0].dataframe()['delta_ln_prob_stdev'].values[1])
    return time, std

def linear_fit(x, b):
    return -0.5*x + b

# after the simulation is complete, perform some tests
class TestFlatHistogramLJ(unittest.TestCase):
    def test(self):
        # analyze grand canonical ensemble average number of particles
        import pandas as pd
        from pyfeasst import macrostate_distribution
        import matplotlib.pyplot as plt
        from scipy.optimize import curve_fit
        show_plot = False
        #show_plot = True
        slopes = list()
        for s in range(3):
            dists = macrostate_distribution.read_appended(file_name='lj_lnpin100s'+str(s)+'.txt', num_states=2)
            #print(dists)
            time, std = dists_to_time_std(dists)
            logt = np.log(time[10:])
            logs = np.log(std[10:])
            popt, pcov = curve_fit(linear_fit, logt, logs)
            slopes.append(popt[0])
            #print(s, popt[0])

            if show_plot:
                label = '2bin'
                if s == 2:
                    label = '1bin'
                plt.plot(np.log(time), np.log(std), label=label)
                plt.plot(logt, linear_fit(logt, popt[0]), label=label + ' slope: ' + str(round(popt[0],2)))

        z12 = np.exp(2*(slopes[2] - np.average(slopes[0:1])))
        print('efficiency of a double macrostate simulation compared to two single macrostate simulations', z12)
        if show_plot:
            plt.xlabel('log(time)', fontsize=16)
            plt.ylabel('log(std_of_average)', fontsize=16)
            #plt.title(r'$z_{12}=$'+str(round(z12, 2)), fontsize=16)
            plt.title('efficiency of a double macrostate simulation\n compared to two single macrostate simulations\n'+r'$z_{12}=$'+str(round(z12, 2)), fontsize=16)
            plt.legend()
            plt.show()

# run the simulation and, if complete, analyze.
def run(sim):
    params['sim'] = sim
    params["seed"] = random.randrange(1e9)
#    if args.task == 0:
    params['window_custom'] = """WindowCustom min0 {min_particles} max {max_particles} overlap 0""".format(**params)
    params['procs_per_sim'] = 1
    if sim >= 2*params['num_1bin_batch']:
        params['window_custom'] = """WindowCustom min0 {min_particles} min1 {max_particles} max {max_particles} overlap 0""".format(**params)
        params['procs_per_sim'] = 1
    file_name = "lj_launchn"+str(params['min_particles'])+"s"+str(sim)+".txt"
    mc_lj(params, file_name=file_name)
    syscode = subprocess.call("~/binning_efficiency/feasst/build/bin/fst < " + file_name + " > lj_launchn"+str(params['min_particles'])+"s"+str(sim)+".log", shell=True, executable='/bin/bash')
#    else:
#        syscode = subprocess.call("~/binning_efficiency/feasst/build/bin/rst lj_checkpointn"+str(params['min_particles'])+"s"+str(sim)+".fst", shell=True, executable='/bin/bash')
    if syscode == 0:
        unittest.main(argv=[''], verbosity=2, exit=False)
    return syscode

if __name__ == "__main__":
    if args.run_type == 0:
        params['num_1bin_batch'] = int(params['procs_per_node']/4)
        with Pool(3*params['num_1bin_batch']) as pool:
            codes = pool.starmap(run, zip(range(0, 3*params['num_1bin_batch'])))
            if np.count_nonzero(codes) > 0:
                sys.exit(1)
    elif args.run_type == 1:
        for min_particles in [200]:
            slurm_queue(min_particles=min_particles)
            subprocess.call("sbatch slurm.txt | awk '{print $4}' >> launch_ids.txt", shell=True, executable='/bin/bash')
    elif args.run_type == 2:
        unittest.main(argv=[''], verbosity=2, exit=False)
    else:
        assert False  # unrecognized run_type
