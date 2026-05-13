"""
Generate a contact table of a protein, by default described by the provided 4lyt.pdb file.
Each atom is modeled as a hard sphere and the contact distance is computed over a number of
orientations.
"""

import copy
import sys
import os.path
import subprocess
import numpy as np
from pyfeasst import fstio
from pyfeasst import coarse_grain_pdb
from launch_1_cg_protein import parse
from launch_1_cg_protein import generate_domain_pairs

def generate_domain_pair(params):
    params['orientation_file'] = 'orientations' + str(params['num_orientations_per_pi']) + '.txt'
    if params['domain1'] != params['domain2']:
        params['particles'] = """particle_type=pt1:{domain1}j{job}.txt,pt2:{domain2}j{job}.txt add_num_pt1_particles=1 add_num_pt2_particles=1""".format(**params)
    else:
        params['particles'] = """particle_type=pt1:{domain1}j{job}.txt add_num_pt1_particles=2""".format(**params)
    return """MonteCarlo
Configuration cubic_side_length=2e2 {particles} cutoff=0 set_cutoff_min_to_sigma=true {vis_extra}
Potential Model=HardSphere VisitModel=VisitModelCell min_length=max_sigma energy_cutoff=1e100
TabulateTwoRigidBody3D proc={sim} num_proc={num_sims} input_orientation_file={orientation_file} num_z=-1 output_table_file={prefix}{sim}_{domain1}_{domain2}_table.txt {contact_xyz_file} contact_xyz_index={contact_xyz_index}
""".format(**params)

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        for domain_pair in params['domain_pairs']:
            #print('domain_pair', domain_pair)
            params['domain1'] = domain_pair[0]
            params['domain2'] = domain_pair[1]
            myfile.write(generate_domain_pair(params))

def post_process(params):
    """ Combining table files, checking length and then launching the next step """
    subprocess.check_call(['sleep', '5']) # without this command, combine doesn't read all tables?
    for domain_pair in params['domain_pairs']:
        #print('domain_pair', domain_pair)
        params['domain1'] = domain_pair[0]
        params['domain2'] = domain_pair[1]
        tabsuffix = '_' + str(domain_pair[0]) + '_' + str(domain_pair[1]) + '_table.txt'
        if not os.path.isfile(params['prefix'] + tabsuffix):
            fstio.combine_tables_two_rigid_body(prefix=params['prefix'], suffix=tabsuffix,
                                                num_procs=params['num_sims'])
            if params['num_orientations_per_pi'] == 2 and params['pH'] == 6 and \
               params['domains'] == '4lyt':
                with open(params['prefix'] + tabsuffix, 'r', encoding="utf-8") as file1:
                    lines = file1.readlines()
                #print(len(lines))
                assert len(lines) == 1125 + 6
                assert lines[0] == 'site_types 1 0\n'
                assert lines[-1] == '-1 160\n'
    print('launching after_1_2_energy.py')
    subprocess.check_call('python after_1_2_energy.py '+fstio.dict_to_argparse(params['original_args']), shell=True, executable='/bin/bash')

def tables_exist(params):
    print('Check if table(s) already exist')
    exists = []
    for domain_pair in params['domain_pairs']:
        #print('domain_pair', domain_pair)
        params['domain1'] = domain_pair[0]
        params['domain2'] = domain_pair[1]
        if os.path.isfile('''{prefix}_{domain1}_{domain2}_table.txt'''.format(**params)):
            exists.append(True)
        else:
            exists.append(False)
    if np.all(exists):
        return True
    return False

if __name__ == '__main__':
    parser = parse()
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    prms = vars(args)
    prms['original_args'] = copy.deepcopy(prms)
    prms['prefix'] = 'contact'
    prms['script'] = __file__
    prms['minutes'] = int(prms['hours_terminate']*60) # minutes allocated on queue
    prms['hours_terminate'] = 0.99*prms['hours_terminate'] - 0.0333 # terminate before queue
    prms['procs_per_sim'] = prms['procs_per_job']
    prms['num_jobs'] = prms['num_jobs_table']
    prms['num_sims'] = prms['num_jobs']*prms['procs_per_job']
    generate_domain_pairs(prms)
    if tables_exist(prms):
        print('Using existing table(s)')
        post_process(prms)
        sys.exit()
    if prms['contact_xyz_file'] != '':
        prms['contact_xyz_file'] = 'contact_xyz_file ' + prms['contact_xyz_file']
        prms['vis_extra'] = 'group=fixed,mobile fixed_particle_type=pt1 mobile_particle_type=pt2'
    else:
        prms['vis_extra'] = ''

    if prms['run_type'] != 2:
        for domain in prms['domains'].split(','):
            print('converting domain:', domain, 'from pdb to pqr')
            # example of how pdb is converted to pqr
            #subprocess.check_call('pdb2pqr30 --titration-state-method=propka --with-ph='+str(prms['pH'])+' --ff=PARSE '+domain+'.pdb '+domain+'.pqr --drop-water > '+domain+'.log 2>&1', shell=True, executable='/bin/bash')
            mol=domain+'j'+str(prms['job'])
            subprocess.check_call('cp '+domain+'.pqr '+mol+'.pqr', shell=True, executable='/bin/bash')
            coarse_grain_pdb.write_fstprt(mol)

    fstio.run_simulations(params=prms,
                          sim_job_dependent_params=None,
                          #sim_job_dependent_params=sim_job_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_job,
                          args=args)
