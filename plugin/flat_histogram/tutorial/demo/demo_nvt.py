"""
In this demonstration of making movies, conduct a canonical ensemble MC simulations of LJ particles.
Usage: python3 demo_nvt.py
When prompted after simulation completes, use VMD to generate a movie as follows
- Extensions->Visualization->Movie Maker
- Renderer -> Internal Tachyon
- Movie Settings -> trajectory, uncheck delete image files
- Format -> MPEG-1 (ffmpeg)
- Set working directory: select current directory of this tutorial
- Click 'Make Movie'
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/particle/lj.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--beta', type=float, default=1./1.5, help='inverse temperature')
PARSER.add_argument('--mu_init', type=float, default=10, help='initial chemical potential')
PARSER.add_argument('--num_particles', type=int, default=370, help='maximum number of particles')
PARSER.add_argument('--cubic_side_length', type=float, default=8,
                    help='cubic periodic boundary length')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e4),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=1e2,
                    help='number of iterations for equilibration')
PARSER.add_argument('--production_iterations', type=int, default=1e3,
                    help='number of iterations for production')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.02, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=0.2, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
PARSER.add_argument('--prefix', type=str, default='lj', help='prefix for all output file names')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--seed', type=int, default=-1,
                    help='Random number generator seed. If -1, assign random seed to each sim.')
PARSER.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
PARSER.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
PARSER.add_argument('--scratch', type=str, default=None,
                    help='Optionally write scheduled job to scratch/logname/jobid.')
PARSER.add_argument('--node', type=int, default=0, help='node ID')
PARSER.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
PARSER.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

# Convert arguments into a parameter dictionary, and add argument-dependent parameters.
ARGS, UNKNOWN_ARGS = PARSER.parse_known_args()
assert len(UNKNOWN_ARGS) == 0, 'An unknown argument was included: '+str(UNKNOWN_ARGS)
PARAMS = vars(ARGS)
PARAMS['script'] = __file__
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.95*PARAMS['hours_terminate'] - 0.05 # terminate FEASST before SLURM
PARAMS['hours_terminate'] *= PARAMS['procs_per_node'] # real time -> cpu time
PARAMS['hours_checkpoint'] *= PARAMS['procs_per_node']
PARAMS['num_sims'] = PARAMS['num_nodes']
PARAMS['procs_per_sim'] = PARAMS['procs_per_node']

def write_feasst_script(params, file_name):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(file_name, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt}
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-4

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Log trials_per_write {trials_per_iteration} file_name {prefix}_eq.txt
Tune
Run until_num_particles {num_particles}
RemoveTrial name TrialAdd
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log

# nvt production
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {production_iterations}
Log trials_per_write {trials_per_iteration} file_name {prefix}.txt
Movie trials_per_write {trials_per_iteration} file_name {prefix}.xyz
Energy trials_per_write {trials_per_iteration} file_name {prefix}_en.txt
Run until_criteria_complete true
""".format(**params))

def post_process(params):
    import subprocess
    print('uncomment the following line to create the simulation image files using vmd')
    #subprocess.check_call('~/software/vmd-1.9.3/build/bin/vmd -e lj.xyz.vmd', shell=True, executable='/bin/bash')
    print('Sometimes, matplotlib requires a command "export QT_QPA_PLATFORM=offscreen" if there are problems with the backend gui')
    log = pd.read_csv(params['prefix']+'.txt')
    print(log['energy'])    
    for frame in range(len(log)):
        fig, ax = plt.subplots()
        # plot the full energy vs trial/frame
        plt.scatter(log['trial'], log['energy'])
        plt.xlabel('trial')
        plt.ylabel('energy')

        # plot a vertical line for the current frame location
        plt.axvline(log['trial'][frame])

        # read the vmd image file and plot as an inset
        # note offset in .xyz frame and log frames. The first xyz frame is discarded (e.g., frame+1)
        im = plt.imread('untitled.'+str("{:05d}".format(frame+1))+'.ppm')
        newax = fig.add_axes([0.17,1,0.5,0.5], anchor='NE', zorder=1)
        newax.imshow(im)
        newax.axis('off')
        
        # save each image
        plt.savefig('demo_nvt'+str("{:05d}".format(frame))+'.png', bbox_inches='tight')

        # clear the figure
        plt.clf()
    # create the video of the animated plot and inset visualization using ffmpeg
    subprocess.check_call('ffmpeg -an -i demo_nvt%05d.png -vcodec mpeg1video -r 24 -qscale:v 10 -vframes 1000 lj_energy.mpg', shell=True, executable='/bin/bash')

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
