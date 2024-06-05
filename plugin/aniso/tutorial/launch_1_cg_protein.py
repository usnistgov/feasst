"""
Generate and simulate a coarse-grained potential of a protein.
This is done in a number of steps.

In the first step, generate an orientation file.
This file is used to determine if the orientation is identical to a previous orientation.
For example, when the spherical polar angle is zero, the orientation is the same for all azimuthal angles.
Orientation files may be included with TabulateTwoRigidBody3D::input_orientation_file.
For each orientation, a value of -1 identifies a unique orientation and a value of zero or more
identifies an earlier orientation that is identical.
These orientation files only need to be generated once per num_orientations_per_pi value.

After this step, the next script after_1_contact is called in post_process.
"""

import argparse
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--num_orientations_per_pi', type=int, default=2,
                    help='maximum number of orientations per 180 degrees')
PARSER.add_argument('--ij', type=int, default=0, help='interaction between different domains.')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--hours_terminate', type=float, default=14*24, help='hours until termination')
PARSER.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
PARSER.add_argument('--procs_per_node', type=int, default=1, help='Number of processors in a node')
PARSER.add_argument('--scratch', type=str, default=None,
                    help='Optionally write scheduled job to scratch/logname/jobid.')
PARSER.add_argument('--queue_flags', type=str, default="", help='extra flags for queue (e.g., for slurm, "-p queue")')
PARSER.add_argument('--node', type=int, default=0, help='node ID')
PARSER.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
PARSER.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

# Convert arguments into a parameter dictionary, and add argument-dependent parameters.
ARGS, UNKNOWN_ARGS = PARSER.parse_known_args()
assert len(UNKNOWN_ARGS) == 0, 'An unknown argument was included: '+str(UNKNOWN_ARGS)
PARAMS = vars(ARGS)
PARAMS['script'] = __file__
PARAMS['prefix'] = 'orientations'
PARAMS['sim_id_file'] = PARAMS['prefix'] + '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.99*PARAMS['hours_terminate'] - 0.0333 # terminate before queue
PARAMS['procs_per_sim'] = PARAMS['procs_per_node']
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']
PARAMS['domain1'] = '/feasst/particle/spce.fstprt'
if PARAMS['ij'] == 1:
    PARAMS['domain2'] = '/feasst/particle/propane.fstprt'
    PARAMS['extra'] = '_ij'
else:
    PARAMS['domain2'] = '/feasst/particle/spce.fstprt'
    PARAMS['extra'] = ''

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
Configuration cubic_side_length 2e2 particle_type0 {domain1} particle_type1 {domain2} \
  add_particles_of_type0 1 add_particles_of_type1 1 \
  group0 fixed fixed_particle_type 0 group1 mobile mobile_particle_type 1
Potential Model ModelEmpty
TabulateTwoRigidBody3D num_orientations_per_pi {num_orientations_per_pi} output_orientation_file {prefix}{num_orientations_per_pi}{extra}.txt
""".format(**params))

def post_process(params):
    nk = params['num_orientations_per_pi']
    expected = (2*nk+1)**2 * (nk+1)**3
    filename = 'orientations' + str(nk)
    if params['ij'] == 1:
        expected = (2*nk+1)**3 * (nk+1)**2
        filename += '_ij'
    filename += '.txt'
    print('expected number of orientations:', expected)
    df = pd.read_csv(filename, skiprows=1, sep='\s+')
    assert expected == len(df.columns)
    print('launching after_1_contact.py')
    subprocess.check_call("""python after_1_contact.py --num_orientations_per_pi {num_orientations_per_pi} --run_type {run_type}""".format(**params), shell=True, executable='/bin/bash')

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
