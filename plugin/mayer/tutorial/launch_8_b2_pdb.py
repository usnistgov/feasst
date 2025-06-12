"""
Compute the second osmotic virial coefficient of rigid all-atom model.
"""

import argparse
import numpy as np
import subprocess
from pyfeasst import fstio
from pyfeasst import physical_constants
from pyfeasst import coarse_grain_pdb

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--temperature', type=float, default=298.15, help='temperature in Kelvin')
    parser.add_argument('--pH', type=float, default=6, help='pH')
    parser.add_argument('--ionic_strength', type=float, default=0.15, help='formulation ionic strength of NaCl in Molar units. If -1, HardSphere only.')
    parser.add_argument('--smoothing_distance', type=float, default=2, help='distance from cutoff to smooth to zero')
    parser.add_argument('--cubic_side_length', type=float, default=1e4, help='cubic side length')
    parser.add_argument('--molecular_weight', type=float, default=14315, help='molecular weight of protein in g/mol')
    parser.add_argument('--num_beta_taylor', type=int, default=10, help='number of Taylor series coefficients')
    parser.add_argument('--reference_sigma', type=float, default=30, help='size of hard sphere on COM of rigid domain')
    parser.add_argument('--trials_per', type=int, default=1e2, help='number of trials per cycle')
    parser.add_argument('--equilibration', type=int, default=1e2, help='number of cycles in equilibration')
    parser.add_argument('--production', type=int, default=1e9, help='number of cycles in production')
    parser.add_argument('--seed', type=int, default=-1,
                        help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--domain1', type=str, default='4lyt', help='fstprt file')
    parser.add_argument('--domain2', type=str, default='4lyt', help='fstprt file')
    parser.add_argument('--run_type', '-r', type=int, default=0,
                        help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--hours_terminate', type=float, default=0.5, help='hours until termination')
    parser.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
    parser.add_argument('--procs_per_node', type=int, default=32, help='Number of nodes in queue')
    parser.add_argument('--scratch', type=str, default=None,
                        help='Optionally write scheduled job to scratch/logname/jobid.')
    parser.add_argument('--queue_flags', type=str, default="", help='extra flags for queue (e.g., for slurm, "-p queue")')
    parser.add_argument('--node', type=int, default=0, help='node ID')
    parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
    parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

    # Convert arguments into a parameter dictionary, and add argument-dependent parameters.
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    params = vars(args)
    params['prefix'] = 'b2pdb'
    params['sim_id_file'] = params['prefix'] + '_sim_ids.txt'
    params['script'] = __file__
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    params['beta'] = 1./(params['temperature']*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
    if params['ionic_strength'] != -1:
        temp_cel = params['temperature'] - 273.15
        params['dielectric_water'] = 87.74 - 0.40008*temp_cel + 9.398e-4*temp_cel**2 - 1.4e-6*temp_cel**3
        eps_0 = physical_constants.VacuumElectricPermittivity().value()
        elem_q = physical_constants.ElementaryCharge().value()
        na = physical_constants.AvogadroConstant().value()
        kb = physical_constants.BoltzmannConstant().value()
        params['kappa'] = np.sqrt(2*(elem_q**2)*params['ionic_strength']*(1e3)*na/(params['dielectric_water']*eps_0*kb*params['temperature']*1e20))
        params['cutoff'] = 5/params['kappa']
        params['cutoff_line'] = """cutoff {cutoff}""".format(**params)
        params['potential'] = """Potential Model=ModelTwoBodyFactory models=HardSphere,LennardJones,DebyeHuckel \
kappa={kappa} dielectric={dielectric_water} smoothing_distance={smoothing_distance} \
energy_cutoff=1e100""".format(**params)
    else:
        params['cutoff_line'] = """set_cutoff_min_to_sigma=true"""
        params['potential'] = """Potential Model=HardSphere energy_cutoff=1e100""".format(**params)

    # convert pdb to pqr to fstprt
    for domain in [params['domain1'], params['domain2']]:
        subprocess.check_call('pdb2pqr30 --titration-state-method=propka --with-ph='+str(params['pH'])+' --ff=PARSE '+domain+'.pdb '+domain+'.pqr --drop-water > '+domain+'.log 2>&1', shell=True, executable='/bin/bash')
        coarse_grain_pdb.write_fstprt(domain)
    if params['queue_id'] == -1:
        with open(params['domain1']+'_with_ref.txt', 'r') as file1:
            params['num_site_types'] = int(file1.readlines()[4].split(' ')[0])
    return params, args

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed={seed}
Configuration cubic_side_length={cubic_side_length} particle_type=domain1:{domain1}_with_ref.txt,domain2:{domain2}_with_ref.txt \
  add_num_domain1_particles=1 add_num_domain2_particles=1 \
  group=mobile mobile_particle_type=domain2 \
  {cutoff_line}
{potential}
RefPotential ref=hs Model=HardSphere sigma=0 cutoff=0 sigma0={reference_sigma} cutoff0={reference_sigma} sigma{num_site_types}={reference_sigma} cutoff{num_site_types}={reference_sigma}
ThermoParams beta={beta}
MayerSampling trials_per_cycle={trials_per} cycles_to_complete={equilibration} num_beta_taylor={num_beta_taylor}
TrialTranslate new_only=true ref=hs tunable_param=35 group=mobile
TrialRotate new_only=true ref=hs tunable_param=40

# tune trial parameters
Let [write]=trials_per_write={trials_per} file_name={prefix}_{domain1}-{domain2}_{sim:03d}
Log [write]_eq.csv
CriteriaWriter [write]_b2_eq.csv
Tune
Run until=complete
Remove name=CriteriaWriter,Log,Tune

# production
Log [write].csv
#Movie [write].xyz
CriteriaWriter [write]_b2.csv
MayerSampling trials_per_cycle={trials_per} cycles_to_complete={production} num_beta_taylor={num_beta_taylor}
Run until=complete
""".format(**params))

def post_process(params):
    from pyfeasst import accumulator
    b2acc = accumulator.Accumulator()
    for p in range(params['procs_per_node']):
        params['sim'] = p
        filename="""{prefix}_{domain1}-{domain2}_{sim:03d}_b2_eq.csv""".format(**params)
        with open(filename, 'r') as file1: lines = file1.readlines()
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
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=None,
                          #sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
