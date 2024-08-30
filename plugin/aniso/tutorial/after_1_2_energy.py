"""
Generate an energy table of a protein, by default described by the provided 4lyt.pdb file.
Each atom is modeled as a hard sphere with an LJ and screened charge interation computed over a number of orientations.
"""

import copy
import sys
import os.path
import subprocess
import numpy as np
from pyfeasst import fstio
from pyfeasst import physical_constants
from launch_1_cg_protein import parse

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
    """ Combine tables, check their length and launch the next step """
    subprocess.check_call(['sleep', '15']) # without this command, combine doesn't read all tables?
    if not os.path.isfile('''{prefix}.txt'''.format(**params)):
        fstio.combine_tables_two_rigid_body(prefix=params['prefix'], suffix='.txt',
                                            num_procs=params['num_sims'])
        if params['num_orientations_per_pi'] == 6 and params['domain1'] == '4lyt' and \
           params['domain2'] == params['domain1']:
            with open("""{prefix}.txt""".format(**params), 'r', encoding="utf-8") as file1:
                lines = file1.readlines()
            #print(len(lines))
            assert len(lines) == 681
            assert lines[0] == 'site_types 1 0\n'
            assert lines[6] == '3.762260e+01 -4.069026e+00 -7.593451e-04\n'
            assert lines[-1] == '-1 160\n'
    print('launching after_1_3_b2.py')
    subprocess.check_call('python after_1_3_b2.py '+fstio.dict_to_argparse(params['original_args']),
                          shell=True, executable='/bin/bash')

if __name__ == '__main__':
    parser = parse()
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    prms = vars(args)
    prms['original_args'] = copy.deepcopy(prms)
    prms['script'] = __file__
    prms['prefix'] = 'energy'
    if os.path.isfile('''{prefix}.txt'''.format(**prms)):
        print('using existing:', '''{prefix}.txt'''.format(**prms))
        post_process(prms)
        sys.exit()
    prms['sim_id_file'] = prms['prefix'] + '_sim_ids.txt'
    prms['minutes'] = int(prms['hours_terminate']*60) # minutes allocated on queue
    prms['hours_terminate'] = 0.99*prms['hours_terminate'] - 0.0333 # terminate before queue
    prms['orientation_file'] = 'orientations' + str(prms['num_orientations_per_pi'])
    if prms['domain1'] != prms['domain2']:
        prms['orientation_file'] += '_ij'
    prms['orientation_file'] += '.txt'
    prms['contact_file'] = 'contact.txt'
    prms['procs_per_sim'] = 1
    prms['num_sims'] = prms['num_nodes']*prms['procs_per_node']
    temp_cel = prms['temperature'] - 273.15
    prms['dielectric_water'] = 87.74 - 0.40008*temp_cel + 9.398e-4*temp_cel**2 - 1.4e-6*temp_cel**3
    eps_0 = physical_constants.VacuumElectricPermittivity().value()
    elem_q = physical_constants.ElementaryCharge().value()
    na = physical_constants.AvogadroConstant().value()
    kb = physical_constants.BoltzmannConstant().value()
    prms['kappa'] = np.sqrt(2*(elem_q**2)*prms['ionic_strength']*(1e3)*na/(prms['dielectric_water']*eps_0*kb*prms['temperature']*1e20))
    prms['cutoff'] = 5/prms['kappa']
    prms['initial_box'] = 4*prms['cutoff'] # initial box adjusted by TabulateTwoRigidBody3D

    fstio.run_simulations(params=prms,
                          sim_node_dependent_params=None,
                          #sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=args)
