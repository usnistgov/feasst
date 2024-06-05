"""
Generate an energy table of a protein, by default described by the provided 4lyt.pdb file.
Each atom is modeled as a hard sphere with an LJ and screened charge interation computed over a number of orientations.
"""

import argparse
import numpy as np
import subprocess
from pyfeasst import fstio
from pyfeasst import physical_constants

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
#PARSER.add_argument('--orientation_file', type=str, default='../orientations/orientations1.txt', help='orientation file')
PARSER.add_argument('--num_orientations_per_pi', type=int, default=1, help='num orientations per 180 degrees')
PARSER.add_argument('--num_z', type=int, default=7, help='num of distances per orientation')
PARSER.add_argument('--gamma', type=float, default=-4, help='stretching exponent for table')
PARSER.add_argument('--temperature', type=float, default=298.15, help='temperature in Kelvin')
PARSER.add_argument('--ionic_strength', type=float, default=0.15, help='formulation ionic strength of NaCl in Molar units')
PARSER.add_argument('--smoothing_distance', type=float, default=2, help='distance from cutoff to smooth to zero')
PARSER.add_argument('--domain1', type=str, default='4lyt', help='fstprt file')
PARSER.add_argument('--domain2', type=str, default='4lyt', help='fstprt file')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--hours_terminate', type=float, default=14*24, help='hours until termination')
PARSER.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
PARSER.add_argument('--procs_per_node', type=int, default=32, help='Number of nodes in queue')
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
PARAMS['prefix'] = 'energy'
PARAMS['sim_id_file'] = PARAMS['prefix'] + '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.99*PARAMS['hours_terminate'] - 0.0333 # terminate before queue
PARAMS['orientation_file'] = 'orientations' + str(PARAMS['num_orientations_per_pi'])
if PARAMS['domain1'] != PARAMS['domain2']:
    PARAMS['orientation_file'] += '_ij'
PARAMS['orientation_file'] += '.txt'
PARAMS['contact_file'] = 'contact.txt'
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']
temp_cel = PARAMS['temperature'] - 273.15
PARAMS['dielectric_water'] = 87.74 - 0.40008*temp_cel + 9.398e-4*temp_cel**2 - 1.4e-6*temp_cel**3
eps_0 = physical_constants.VacuumElectricPermittivity().value()
elem_q = physical_constants.ElementaryCharge().value()
na = physical_constants.AvogadroConstant().value()
kb = physical_constants.BoltzmannConstant().value()
PARAMS['kappa'] = np.sqrt(2*(elem_q**2)*PARAMS['ionic_strength']*(1e3)*na/(PARAMS['dielectric_water']*eps_0*kb*PARAMS['temperature']*1e20))
PARAMS['cutoff'] = 5/PARAMS['kappa']
PARAMS['initial_box'] = 4*PARAMS['cutoff'] # assume initial box to fit cutoff, but TabulateTwoRigidBody3D will adjust it

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
Configuration cubic_side_length {initial_box} particle_type0 {domain1}.fstprt particle_type1 {domain2}.fstprt \
  add_particles_of_type0 1 add_particles_of_type1 1 cutoff {cutoff} \
  group0 fixed fixed_particle_type 0 group1 mobile mobile_particle_type 1
Potential Model ModelTwoBodyFactory \
  model0 HardSphere \
  model1 LennardJones \
  model2 DebyeHuckel kappa {kappa} dielectric {dielectric_water} smoothing_distance {smoothing_distance} \
  VisitModel VisitModelCell min_length max_cutoff energy_cutoff 1e100
TabulateTwoRigidBody3D proc {sim} num_proc {num_sims} input_orientation_file {orientation_file} num_z {num_z} smoothing_distance {smoothing_distance} input_table_file {contact_file} output_table_file {prefix}{sim}.txt gamma {gamma}
""".format(**params))

def post_process(params):
    subprocess.check_call(['sleep', '5'])
    fstio.combine_tables_two_rigid_body(prefix=params['prefix'], suffix='.txt', num_procs=params['num_sims'])
    if params['num_orientations_per_pi'] == 6 and params['domain1'] == '4lyt' and params['domain2'] == params['domain1']:
        with open("""{prefix}.txt""".format(**params), 'r') as file1:
            lines = file1.readlines()
        #print(len(lines))
        assert len(lines) == 681
        assert lines[0] == 'site_types 1 0\n'
        assert lines[6] == '3.762260e+01 -4.069026e+00 -7.593451e-04\n'
        assert lines[-1] == '-1 160\n'
    print('launching after_1_b2.py')
    subprocess.check_call("""python after_1_b2.py --run_type {run_type}""".format(**params), shell=True, executable='/bin/bash')

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          #sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
