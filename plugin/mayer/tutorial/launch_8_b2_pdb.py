"""
Compute the second osmotic virial coefficient of rigid all-atom model.
"""

import argparse
import numpy as np
import subprocess
from pyfeasst import fstio
from pyfeasst import physical_constants
from pyfeasst import coarse_grain_pdb

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--temperature', type=float, default=298.15, help='temperature in Kelvin')
PARSER.add_argument('--pH', type=float, default=6, help='pH')
PARSER.add_argument('--ionic_strength', type=float, default=0.15, help='formulation ionic strength of NaCl in Molar units. If -1, HardSphere only.')
PARSER.add_argument('--smoothing_distance', type=float, default=2, help='distance from cutoff to smooth to zero')
PARSER.add_argument('--cubic_side_length', type=float, default=1e4, help='cubic side length')
PARSER.add_argument('--molecular_weight', type=float, default=14315, help='molecular weight of protein in g/mol')
PARSER.add_argument('--num_beta_taylor', type=int, default=10, help='number of Taylor series coefficients')
PARSER.add_argument('--reference_sigma', type=float, default=30, help='size of hard sphere on COM of rigid domain')
PARSER.add_argument('--trials_per', type=int, default=1e2, help='number of trials per iteration')
PARSER.add_argument('--equilibration', type=int, default=1e2, help='number of iterations in equilibration')
PARSER.add_argument('--production', type=int, default=1e9, help='number of iterations in production')
PARSER.add_argument('--seed', type=int, default=-1,
                    help='Random number generator seed. If -1, assign random seed to each sim.')
PARSER.add_argument('--domain1', type=str, default='4lyt', help='fstprt file')
PARSER.add_argument('--domain2', type=str, default='4lyt', help='fstprt file')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--hours_terminate', type=float, default=0.5, help='hours until termination')
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
PARAMS['prefix'] = 'b2pdb'
PARAMS['sim_id_file'] = PARAMS['prefix'] + '_sim_ids.txt'
PARAMS['script'] = __file__
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.99*PARAMS['hours_terminate'] - 0.0333 # terminate before queue
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']
PARAMS['beta'] = 1./(PARAMS['temperature']*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
if PARAMS['ionic_strength'] != -1:
    temp_cel = PARAMS['temperature'] - 273.15
    PARAMS['dielectric_water'] = 87.74 - 0.40008*temp_cel + 9.398e-4*temp_cel**2 - 1.4e-6*temp_cel**3
    eps_0 = physical_constants.VacuumElectricPermittivity().value()
    elem_q = physical_constants.ElementaryCharge().value()
    na = physical_constants.AvogadroConstant().value()
    kb = physical_constants.BoltzmannConstant().value()
    PARAMS['kappa'] = np.sqrt(2*(elem_q**2)*PARAMS['ionic_strength']*(1e3)*na/(PARAMS['dielectric_water']*eps_0*kb*PARAMS['temperature']*1e20))
    PARAMS['cutoff'] = 5/PARAMS['kappa']
    PARAMS['cutoff_line'] = """cutoff {cutoff}""".format(**PARAMS)
    PARAMS['potential'] = """Potential Model ModelTwoBodyFactory \
model0 HardSphere \
model1 LennardJones \
model2 DebyeHuckel kappa {kappa} dielectric {dielectric_water} smoothing_distance {smoothing_distance} \
energy_cutoff 1e100""".format(**PARAMS)
else:
    PARAMS['cutoff_line'] = """set_cutoff_min_to_sigma true"""
    PARAMS['potential'] = """Potential Model HardSphere energy_cutoff 1e100""".format(**PARAMS)

# convert pdb to pqr to fstprt
for domain in [PARAMS['domain1'], PARAMS['domain2']]:
    subprocess.check_call('pdb2pqr30 --titration-state-method=propka --with-ph='+str(PARAMS['pH'])+' --ff=PARSE '+domain+'.pdb '+domain+'.pqr --drop-water > '+domain+'.log 2>&1', shell=True, executable='/bin/bash')
    coarse_grain_pdb.write_fstprt(domain)
if PARAMS['queue_id'] == -1:
    with open(PARAMS['domain1']+'_with_ref.fstprt', 'r') as file1:
        PARAMS['num_site_types'] = int(file1.readlines()[4].split(' ')[0])

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {domain1}_with_ref.fstprt particle_type1 {domain2}_with_ref.fstprt \
  add_particles_of_type0 1 add_particles_of_type1 1 \
  group0 mobile mobile_particle_type 1 \
  {cutoff_line}
{potential}
RefPotential Model HardSphere sigma 0 cutoff 0 sigma0 {reference_sigma} cutoff0 {reference_sigma} sigma{num_site_types} {reference_sigma} cutoff{num_site_types} {reference_sigma}
ThermoParams beta {beta}
MayerSampling num_trials_per_iteration {trials_per} num_iterations_to_complete {equilibration} num_beta_taylor {num_beta_taylor}
TrialTranslate new_only true reference_index 0 tunable_param 35 group mobile
TrialRotate new_only true reference_index 0 tunable_param 40

# tune trial parameters
CriteriaWriter trials_per_write {trials_per} file_name {prefix}_{domain1}-{domain2}_{sim}_b2_eq.txt
Log trials_per_write {trials_per} file_name {prefix}_{domain1}-{domain2}_{sim}_eq.txt
#Movie trials_per_write {trials_per} file_name {prefix}_{domain1}-{domain2}_{sim}_eq.xyz
#CPUTime trials_per_write {trials_per} file_name {prefix}_{domain1}-{domain2}_{sim}_cpu.txt append true
Tune
#Run num_trials {equilibration}
Run until_criteria_complete true
RemoveAnalyze name CriteriaWriter
RemoveAnalyze name Log
#RemoveAnalyze name Movie
RemoveModify name Tune

# production
CriteriaWriter trials_per_write {trials_per} file_name {prefix}_{domain1}-{domain2}_{sim}_b2.txt
Log trials_per_write {trials_per} file_name {prefix}_{domain1}-{domain2}_{sim}.txt
#Movie trials_per_write {trials_per} file_name {prefix}_{domain1}-{domain2}_{sim}.xyz
MayerSampling num_trials_per_iteration {trials_per} num_iterations_to_complete {production} num_beta_taylor {num_beta_taylor}
Run until_criteria_complete true
""".format(**params))

def post_process(params):
    from pyfeasst import accumulator
    b2acc = accumulator.Accumulator()
    for p in range(params['procs_per_node']):
        with open(params['prefix']+'_'+params['domain1']+'-'+params['domain2']+'_'+str(p)+"_b2_eq.txt", 'r') as file1:
            lines = file1.readlines()
        print('p', p, 'lines[0]', lines[0])
        exec('iprm=' + lines[0], globals())
        #print(iprm)
        ref = (2*np.pi/3)*params['reference_sigma']**3
        #print('ref(A^3)', ref)
        b2=iprm['second_virial_ratio']*ref
        #print('b2(A^3)', b2)
        mw = params['molecular_weight']/1e3 # kDa
        b2_molmlg2 = b2*1e-26*physical_constants.AvogadroConstant().value()/mw/mw
        #print(b2_molmlg2, '10^-4 mol*ml*g^-2')
        b2acc.add(b2_molmlg2)
    print('b2', b2acc.mean(), b2acc.stdev()/np.sqrt(params['procs_per_node']))

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          #sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
