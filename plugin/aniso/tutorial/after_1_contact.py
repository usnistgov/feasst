"""
Generate a contact table of a protein, by default described by the provided 4lyt.pdb file.
Each atom is modeled as a hard sphere and the contact distance is computed over a number of orientations.
"""

import argparse
import subprocess
from pyfeasst import fstio
from pyfeasst import coarse_grain_pdb

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--num_orientations_per_pi', type=int, default=2, help='num orientations per 180 degrees')
PARSER.add_argument('--pH', type=float, default=6, help='pH')
PARSER.add_argument('--domain1', type=str, default='4lyt', help='fstprt file')
PARSER.add_argument('--domain2', type=str, default='4lyt', help='fstprt file')
PARSER.add_argument('--contact_xyz_file', type=str, default='', help='If not empty, print contact configuration for each orientation')
PARSER.add_argument('--contact_xyz_index', type=int, default=-1, help='If not -1, print contact configuration for one orientation')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--hours_terminate', type=float, default=14*24, help='hours until termination')
PARSER.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
PARSER.add_argument('--procs_per_node', type=int, default=16, help='Number of nodes in queue')
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
PARAMS['prefix'] = 'contact'
PARAMS['script'] = __file__
PARAMS['sim_id_file'] = PARAMS['prefix'] + '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.99*PARAMS['hours_terminate'] - 0.0333 # terminate before queue
PARAMS['orientation_file'] = 'orientations' + str(PARAMS['num_orientations_per_pi'])
if PARAMS['domain1'] != PARAMS['domain2']:
    PARAMS['orientation_file'] += '_ij'
PARAMS['orientation_file'] += '.txt'
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']
if PARAMS['contact_xyz_file'] != '':
    PARAMS['contact_xyz_file'] = 'contact_xyz_file ' + PARAMS['contact_xyz_file']
    PARAMS['vis_extra'] = 'group0 fixed fixed_particle_type 0 group1 mobile mobile_particle_type 1'
else:
    PARAMS['vis_extra'] = ''

# convert pdb to pqr to fstprt
if PARAMS['run_type'] != 2:
    for domain in [PARAMS['domain1'], PARAMS['domain2']]:
        subprocess.check_call('pdb2pqr30 --titration-state-method=propka --with-ph='+str(PARAMS['pH'])+' --ff=PARSE '+domain+'.pdb '+domain+'.pqr --drop-water > '+domain+'.log 2>&1', shell=True, executable='/bin/bash')
        coarse_grain_pdb.write_fstprt(domain)

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
Configuration cubic_side_length 2e2 particle_type0 {domain1}.fstprt particle_type1 {domain2}.fstprt \
  add_particles_of_type0 1 add_particles_of_type1 1 set_cutoff_min_to_sigma true {vis_extra}
Potential Model HardSphere VisitModel VisitModelCell min_length max_sigma energy_cutoff 1e100
TabulateTwoRigidBody3D proc {sim} num_proc {num_sims} input_orientation_file {orientation_file} num_z -1 output_table_file {prefix}{sim}.txt {contact_xyz_file} contact_xyz_index {contact_xyz_index}
""".format(**params))

def post_process(params):
    subprocess.check_call(['sleep', '5'])
    fstio.combine_tables_two_rigid_body(prefix=params['prefix'], suffix='.txt', num_procs=params['num_sims'])
    if params['num_orientations_per_pi'] == 6 and params['pH'] == 6 and params['domain1'] == '4lyt' and params['domain2'] == params['domain1']:
        with open("""{prefix}.txt""".format(**params), 'r') as file1:
            lines = file1.readlines()
        #print(len(lines))
        assert len(lines) == 681
        assert lines[0] == 'site_types 1 0\n'
        assert lines[-1] == '-1 160\n'
    print('launching after_1_energy.py')
    subprocess.check_call("""python after_1_energy.py --num_orientations_per_pi {num_orientations_per_pi} --run_type {run_type}""".format(**params), shell=True, executable='/bin/bash')

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          #sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
