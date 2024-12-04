"""
Compute the second virial coefficient using coarse-grained models.
"""

import subprocess
import numpy as np
import pandas as pd
from pyfeasst import accumulator
from pyfeasst import fstio
from pyfeasst import physical_constants
from launch_1_cg_protein import parse

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt} add_particles_of_type0 2 group0 mobile mobile_particle_index 1
Potential Model TwoBodyTable VisitModelInner VisitModelInnerTable table_file {table_file} ignore_energy {ignore_energy}
RefPotential Model HardSphere sigma 0 cutoff 0 sigma0 {reference_sigma} cutoff0 {reference_sigma}
ThermoParams beta {beta}
MayerSampling num_beta_taylor {num_beta_taylor} trials_per_cycle {trials_per} cycles_to_complete {equilibration}
TrialTranslate new_only true reference_index 0 tunable_param 35 group mobile
TrialRotate new_only true reference_index 0 tunable_param 40

# tune trial parameters
CriteriaWriter trials_per_write {trials_per} output_file {prefix}_{sim}_b2_eq.txt
Log trials_per_write {trials_per} output_file {prefix}_{sim}_eq.txt
Tune
Run until complete
Remove name0 CriteriaWriter name1 Log name2 Tune

# production
CriteriaWriter trials_per_write {trials_per} output_file {prefix}_{sim}_b2.txt
Log trials_per_write {trials_per} output_file {prefix}_{sim}.txt
#Movie trials_per_write {trials_per} output_file {prefix}_{sim}.xyz
MayerSampling num_beta_taylor {num_beta_taylor} trials_per_cycle {trials_per} cycles_to_complete {production}
Run until complete
""".format(**params))

def post_process(params):
    """ Compute statistics and output B2 value in appropriate units. """
    subprocess.check_call(['sleep', '5']) # without this command, combine doesn't read all tables?
    b2acc = accumulator.Accumulator()
    ref = (2*np.pi/3)*params['reference_sigma']**3
    mw = params['molecular_weight']/1e3 # kDa
    ref *= 1e-26*physical_constants.AvogadroConstant().value()/mw/mw # A3 to 10^4 molml/g2
    for p in range(params['procs_per_node']):
        with open(params['prefix']+'_'+str(p)+"_b2.txt", 'r', encoding="utf-8") as file1:
#        with open(params['prefix']+'_'+str(p)+"_b2_eq.txt", 'r') as file1:
            lines = file1.readlines()
        if len(lines) != 0:
            #print('p', p, 'lines[0]', lines[0])
            exec('iprm=' + lines[0], globals())
            #print(iprm)
            #print('ref(A^3)', ref)
            b2_molmlg2=iprm['second_virial_ratio']*ref
            #print('b2(A^3)', b2)
            #print(b2_molmlg2, '10^-4 mol*ml*g^-2')
            b2acc.add(b2_molmlg2)
            if p == 0:
                df = pd.DataFrame(data={p: iprm['beta_taylor']})
            else:
                df[p] = iprm['beta_taylor']
    print('b2 (mol*ml/g2)', b2acc.mean(), b2acc.stdev()/np.sqrt(params['procs_per_node']))
    if params['molecular_weight'] == 14315.03534:
        assert np.abs(b2acc.mean() - 1.72) < 0.1

if __name__ == '__main__':
    parser = parse()
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    prms = vars(args)
    prms['prefix'] = 'b2'
    prms['sim_id_file'] = prms['prefix'] + '_sim_ids.txt'
    prms['script'] = __file__
    prms['minutes'] = int(prms['hours_terminate']*60) # minutes allocated on queue
    prms['hours_terminate'] = 0.99*prms['hours_terminate'] - 0.0333 # terminate before queue
    prms['procs_per_sim'] = 1
    prms['num_sims'] = prms['num_nodes']*prms['procs_per_node']
    prms['beta'] = 1./(prms['temperature']*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ

    fstio.run_simulations(params=prms,
                          sim_node_dependent_params=None,
                          #sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=args)
