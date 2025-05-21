import sys
import subprocess
import argparse
import random

# define parameters of a pure component NVT MC Lennard-Jones simulation
params = {
    "cubic_side_length": 5, "fstprt": "/feasst/particle/lj.txt", "beta": 1/1.2,
    #"cubic_side_length": 6, "fstprt": "/feasst/particle/lj.txt", "beta": 1/1.5,
    "max_particles": 90, "min_particles": 0, "min_sweeps": 1, "mu": -2.352321,
    #"max_particles": 156, "min_particles": 0, "min_sweeps": 1e4, "mu": -2.352321,
    "trials_per": int(5e3), "hours_per_adjust": 0.01, "hours_per_checkpoint": 1, "seed": random.randrange(int(1e9)), "num_hours": 5*24,
    #"trials_per": int(5e3), "hours_per_adjust": 0.01, "hours_per_checkpoint": 1, "seed": random.randrange(int(1e9)), "num_hours": 5*24,
    "equilibration": 1e6, "num_nodes": 1, "procs_per_node": 1, "script": __file__, "dccb_cut": 2**(1./6.),
    "min_window_size": 5}
params["num_minutes"] = round(params["num_hours"]*60)
params["hours_per_adjust"] = params["hours_per_adjust"]*params["procs_per_node"]
params["hours_per_checkpoint"] = params["hours_per_checkpoint"]*params["procs_per_node"]
params["num_hours_terminate"] = 0.95*params["num_hours"]*params["procs_per_node"]
params["mu_init"] = params["mu"] + 1

# write fst script to run a single simulation
def mc_lj(params=params, file_name="launch.txt"):
    with open(file_name, "w") as myfile: myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_per_adjust} ln_prob_file lj_lnpi.txt bounds_file lj_bounds.txt num_adjust_per_write 10 min_window_size {min_window_size} ln_prob_file_append true
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 2.25 min_size {min_window_size}
Checkpoint file_name lj_checkpoint.fst num_hours {hours_per_checkpoint} num_hours_terminate {num_hours_terminate}

# begin description of each MC clone
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt} cutoff 2.5
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
RefPotential Model LennardJones VisitModel VisitModelCell min_length {dccb_cut}
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
Tune
CheckEnergy trials_per_update {trials_per} tolerance 1e-8
#Checkpoint file_name lj_checkpoint[sim_index].fst num_hours {hours_per_checkpoint}

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Run until_num_particles [soft_macro_min]
Remove name TrialAdd
Run num_trials {equilibration}
Remove name Tune

# gcmc tm production
ThermoParams beta {beta} chemical_potential {mu}
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias TransitionMatrix min_sweeps {min_sweeps} new_sweep 1
#Bias WLTM min_sweeps {min_sweeps} new_sweep 1 min_flatness 25 collect_flatness 20 min_collect_sweeps 20
TrialTransfer weight 2 particle_type 0 reference_index 0 num_steps 4
Movie trials_per_write {trials_per} file_name lj[sim_index].xyz
Log trials_per_write {trials_per} file_name lj[sim_index].txt
Tune trials_per_write {trials_per} file_name lj_tune[sim_index].txt multistate true stop_after_cycle 1
Energy trials_per_write {trials_per} file_name lj_en[sim_index].txt multistate true start_after_cycle 1
CriteriaUpdater trials_per_update {trials_per}
CriteriaWriter trials_per_write {trials_per} file_name lj_crit[sim_index].txt append true
#Run until complete
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
parser.add_argument('--run_type', '-r', type=int, default=1, help="0: submit batch to scheduler, 1: run batch on host")
parser.add_argument('--task', type=int, default=0, help="input by slurm scheduler. If >0, restart from checkpoint.")
args = parser.parse_args()

# after the simulation is complete, perform some tests
def test():
    # analyze grand canonical ensemble average number of particles
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from pyfeasst import macrostate_distribution
#        dists = macrostate_distribution.read_appended('lj_crit00.txt', num_states=params['max_particles']+1)
#        print(dists[0][0].ln_prob())
    with open('lj_crit00.txt', 'r') as file1:
        lines = file1.readlines()
    lines_per_frame = params['max_particles']+4
    print('lines_per_frame', lines_per_frame)
    frames = int(round(len(lines)/lines_per_frame))
    print('len', len(lines))
    print('frames', frames)
    num = pd.read_csv('lj00.txt')
    skip = 1
    for frame in range(1, frames, skip):
        print('frame', frame)
        start = lines_per_frame*frame
        end = lines_per_frame*(frame+1)
        tmpfile = 'ljtmp'+str(frame)
        with open(tmpfile, 'w') as file1:
            for line in lines[start:end]:
                file1.write(line)
        lnpi = macrostate_distribution.MacrostateDistribution(file_name=tmpfile)
        lnpi.reweight(-1.0743260638585548, inplace=True)
        #if frame == frames - 1:
            #lnpi.reweight(-1, inplace=True)
            #print(lnpi.equilibrium())
        fac = (np.pi/6)/params['cubic_side_length']**3
        volfrac = lnpi.macrostates()*fac
        prob = np.exp(lnpi.ln_prob())
        fig, ax = plt.subplots()
        ax.scatter(volfrac, prob)
        #print('num', num, num['p0'])
        ax.scatter(num['num_particles_of_type0'][frame]*fac, 0.001, color='red')
        im = plt.imread('tmp/untitled.'+str("{:05d}".format(frame))+'.ppm')
        plt.xlabel('Volume Fraction', fontsize=20)
        plt.ylabel('Probability', fontsize=20)
        #plt.ylim([-12, -2])
        plt.ylim([0, 0.065])
        newax = fig.add_axes([0.17,0.5,0.5,0.5], anchor='NE', zorder=1)
        newax.imshow(im)
        newax.axis('off')
        plt.savefig('lj'+str("{:04d}".format(frame))+'.png', bbox_inches='tight')
        #plt.show()
        plt.close()
        plt.clf()
        # ffmpeg -i lj%04d.png -c:v libx264 -crf 23 -profile:v baseline -level 3.0 -pix_fmt yuv420p -c:a aac -ac 2 -b:a 128k -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -movflags faststart output.mp4
        # ffmpeg -an -i lj%04d.png -vcodec mpeg1video -r 24 -qscale:v 10 -vframes 500 lnpi.mpg
        # ffmpeg -an -i lj%04d.png -vcodec libx264 -pix_fmt yuv420p -r 24 -qscale:v 10 -vframes 500 lnpi.mpg
        # not this#ffmpeg -an -i /home/user/feasst/plugin/flat_histogram/tutorial/demo/tmp/untitled.%05d.ppm -vcodec mpeg1video -r 24 -vframes 2967 /home/user/feasst/plugin/flat_histogram/tutorial/demo/tmp/untitled.mpg
        # convert old: ffmpeg -i gca.mpg -c:v libx264 -crf 23 -profile:v baseline -level 3.0 -pix_fmt yuv420p -c:a aac -ac 2 -b:a 128k -movflags faststart -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" gca2.mp4

# run the simulation and, if complete, analyze.
def run():
    if args.task == 0:
        file_name = "lj_launch.txt"
        mc_lj(params, file_name=file_name)
        syscode = subprocess.call("../../../../build/bin/fst < " + file_name + " > lj_launch.log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../../../build/bin/rst lj_checkpoint.fst", shell=True, executable='/bin/bash')
    if syscode == 0:
        test()
    return syscode

if __name__ == "__main__":
    if args.run_type == 1:
        syscode = run()
        if syscode != 0:
            sys.exit(1)
    if args.run_type == 2:
        test()
    else:
        assert False  # unrecognized run_type
