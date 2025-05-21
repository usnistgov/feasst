"""
Generate and simulate a coarse-grained potential of a protein.
This is done in a number of steps.

In the first step, generate an orientation file.
This file is used to determine if the orientation is identical to a previous orientation.
For example, when the spherical polar angle is zero, the orientation is the same for all azimuthal
angles.
Orientation files may be included with TabulateTwoRigidBody3D::input_orientation_file.
For each orientation, a value of -1 identifies a unique orientation and a value of zero or more
identifies an earlier orientation that is identical.
These orientation files only need to be generated once per num_orientations_per_pi value.

After this step, the next script after_1_contact is called in post_process.
"""

import copy
import os.path
import argparse
import subprocess
import pandas as pd
from pyfeasst import fstio

def parse():
    """ Parse arguments for this and following scripts """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--num_orientations_per_pi', type=int, default=2,
                        help='maximum number of orientations per 180 degrees')
    parser.add_argument('--run_type', '-r', type=int, default=0,
                        help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--hours_terminate', type=float, default=14*24,
                        help='hours until termination')
    parser.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
    parser.add_argument('--procs_per_node', type=int, default=16, help='Number processors in a node')
    parser.add_argument('--scratch', type=str, default=None,
                        help='Optionally write scheduled job to scratch/logname/jobid.')
    parser.add_argument('--queue_flags', type=str, default="",
                        help='extra flags for queue (e.g., for slurm, "-p queue")')
    parser.add_argument('--node', type=int, default=0, help='node ID')
    parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
    parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

    # additional in contact step (or more)
    parser.add_argument('--pH', type=float, default=6, help='pH')
    parser.add_argument('--domain1', type=str, default='4lyt', help='fstprt file')
    parser.add_argument('--domain2', type=str, default='4lyt', help='fstprt file')
    parser.add_argument('--contact_xyz_file', type=str, default='',
                        help='If not empty, print contact configuration for each orientation')
    parser.add_argument('--contact_xyz_index', type=int, default=-1,
                        help='If not -1, print contact configuration for one orientation')

    # additional in energy step (or more)
    parser.add_argument('--num_z', type=int, default=7, help='num of distances per orientation')
    parser.add_argument('--gamma', type=float, default=-4, help='stretching exponent for table')
    parser.add_argument('--temperature', type=float, default=298.15, help='temperature in Kelvin')
    parser.add_argument('--ionic_strength', type=float, default=0.15, help='formulation ionic strength of NaCl in Molar units')
    parser.add_argument('--smoothing_distance', type=float, default=2, help='distance from cutoff to smooth to zero')

    # additional in b2 step
    parser.add_argument('--table_file', type=str, default='energy.txt', help='file describe cg table potential.')
    parser.add_argument('--cubic_side_length', type=float, default=1e4, help='cubic side length')
    parser.add_argument('--molecular_weight', type=float, default=14315.03534, help='molecular weight of protein in g/mol')
    parser.add_argument('--reference_sigma', type=float, default=30, help='size of hard sphere on COM of rigid domain')
    parser.add_argument('--ignore_energy', type=str, default='false', help='true if interaction is excluded volume only.')
    parser.add_argument('--ignore_intra_energy', type=str, default='false', help='true if intra interaction is excluded volume only.')
    parser.add_argument('--num_beta_taylor', type=int, default=10, help='number of Tayler series derivatives')
    parser.add_argument('--show_plot', type=int, default=0, help='show extrapolation plot if != 0')
    parser.add_argument('--trials_per', type=int, default=int(1e5), help='number of trials per cycle')
    parser.add_argument('--equilibration', type=int, default=int(2e1), help='number of cycles in equilibration')
    parser.add_argument('--production', type=int, default=int(2e1), help='number of cycles in production')
    parser.add_argument('--seed', type=int, default=-1,
                        help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--fstprt', type=str, default='/feasst/plugin/aniso/particle/aniso_tabular.txt', help='fstprt file')
    return parser

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
Configuration cubic_side_length 2e2 particle_type0 /feasst/particle/spce.txt particle_type1 /feasst/particle/spce.txt \
  add_particles_of_type0 1 add_particles_of_type1 1 \
  group0 fixed fixed_particle_type 0 group1 mobile mobile_particle_type 1
Potential Model ModelEmpty
TabulateTwoRigidBody3D num_orientations_per_pi {num_orientations_per_pi} output_orientation_file {prefix}{num_orientations_per_pi}.txt

MonteCarlo
Configuration cubic_side_length 2e2 particle_type0 /feasst/particle/spce.txt particle_type1 /feasst/particle/propane.txt \
  add_particles_of_type0 1 add_particles_of_type1 1 \
  group0 fixed fixed_particle_type 0 group1 mobile mobile_particle_type 1
Potential Model ModelEmpty
TabulateTwoRigidBody3D num_orientations_per_pi {num_orientations_per_pi} output_orientation_file {prefix}{num_orientations_per_pi}_ij.txt
""".format(**params))

def post_process(params):
    """ Check the final file length and then launch the next step. """
    nk = params['num_orientations_per_pi']
    for ij in [True, False]:
        extra = ''
        expected = (2*nk+1)**2 * (nk+1)**3
        if ij:
            expected = (2*nk+1)**3 * (nk+1)**2
            extra = '_ij'
        #print('expected number of orientations:', expected)
        df = pd.read_csv('''{prefix}{num_orientations_per_pi}'''.format(**params)+extra+'.txt', skiprows=1, sep=r'\s+')
        assert expected == len(df.columns)
    print('launching after_1_1_contact.py')
    subprocess.check_call('python after_1_1_contact.py '+fstio.dict_to_argparse(params['original_args']),
                          shell=True, executable='/bin/bash')

if __name__ == '__main__':
    prsr = parse()
    args, unknown_args = prsr.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    prms = vars(args)
    prms['original_args'] = copy.deepcopy(prms)
    prms['script'] = __file__
    prms['prefix'] = 'orientations'
    prms['sim_id_file'] = prms['prefix'] + '_sim_ids.txt'
    prms['minutes'] = int(prms['hours_terminate']*60) # minutes allocated on queue
    prms['hours_terminate'] = 0.99*prms['hours_terminate'] - 0.0333 # terminate before queue
    prms['procs_per_sim'] = prms['procs_per_node']
    prms['num_nodes'] = 1
    prms['num_sims'] = prms['num_nodes']
    orient_file_prefix = '''{prefix}{num_orientations_per_pi}'''.format(**prms)
    if os.path.isfile(orient_file_prefix+'.txt') and \
       os.path.isfile(orient_file_prefix+'_ij.txt'):
        print('using existing:', orient_file_prefix)
        post_process(params=prms)
    else:
        fstio.run_simulations(params=prms,
                              sim_node_dependent_params=None,
                              write_feasst_script=write_feasst_script,
                              post_process=post_process,
                              queue_function=fstio.slurm_single_node,
                              args=args)
