"""
This tutorial is similar to tutorial 5, but for this low temperature simulation,
we will split the simulation into two different nodes.
The first node will have less particles but a higher number of sweeps required.
The second node will have dccb.
This tutorial uses flat histogram LJ tutorial for most of the simulation setup.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import physical_constants
from pyfeasst import macrostate_distribution
from pyfeasst import multistate_accumulator
import launch_04_lj_tm_parallel

def parse(temperature=300):
    """ Parse arguments from command line or change their default values. """
    params, args = launch_04_lj_tm_parallel.parse(
          fstprt='/feasst/particle/spce.fstprt',
          beta=1./(temperature*physical_constants.MolarGasConstant().value()/1e3), # mol/kJ
          beta_mu=-15.24,
          mu_init=-7,
          min_sweeps=5,
          cubic_side_length=20,
          max_particles=296,
          window_alpha=1.5,
          hours_checkpoint=1,
          hours_terminate=5*24,
          num_nodes=2,
          collect_flatness=18,
          min_flatness=22)
    params['script'] = __file__
    params['prefix'] = 'spce_lowt'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['dccb_cut'] = 0.75*3.165
    params['ewald_alpha'] = 5.6/params['cubic_side_length']
    params['num_particles_first_node'] = 180

    # write TrialGrowFile
    with open(params['prefix']+'_grow.txt', 'w') as f:
        f.write("""TrialGrowFile

particle_type 0 weight 2 transfer true site 0 num_steps 10 reference_index 0
bond true mobile_site 1 anchor_site 0 reference_index 0
angle true mobile_site 2 anchor_site 0 anchor_site2 1 reference_index 0
""")

    return params, args

def sim_node_dependent_params(params):
    """ Define parameters that are dependent on the sim or node. """
    if params['node'] == 0:
        params['min_particles'] = 0
        params['max_particles'] = params['num_particles_first_node']
        params['muvt_trials'] = 'TrialTransfer weight 2 particle_type 0'
        params['ref_potential'] = ''
        params['min_sweeps'] = 10
        params['window_alpha'] = 1.1
        params['min_window_size'] = 5
    elif params['node'] == 1:
        params['min_particles'] = params['num_particles_first_node']
        params["muvt_trials"]="""TrialGrowFile grow_file {prefix}_grow.txt""".format(**params)
        params["ref_potential"]="""RefPotential Model HardSphere group oxygen cutoff {dccb_cut} VisitModel VisitModelCell min_length {dccb_cut} cell_group oxygen""".format(**params)
        params['min_sweeps'] = 1
        params['window_alpha'] = 1.25
        params['min_window_size'] = 3
    params['system'] = """Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt} group0 oxygen oxygen_site_type 0
Potential VisitModel Ewald alpha {ewald_alpha} kmax_squared 38
Potential Model ModelTwoBodyFactory model0 LennardJones model1 ChargeScreened erfc_table_size 2e4 VisitModel VisitModelCutoffOuter
{ref_potential}
Potential Model ChargeScreenedIntra VisitModel VisitModelBond
Potential Model ChargeSelf
Potential VisitModel LongRangeCorrections""".format(**params)
    params['nvt_trials'] = """TrialTranslate weight 0.5 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 0.5 particle_type 0 tunable_param 0.5 tunable_target_acceptance 0.25"""
    params['sim_start'] = params['node']*params['procs_per_node']
    params['sim_end'] = params['sim_start'] + params['procs_per_node'] - 1
    params['windows'] = macrostate_distribution.window_exponential(
        alpha=params['window_alpha'], minimums=[params['min_particles']], maximum=params['max_particles'],
        number=params['procs_per_node'], overlap=1, min_size=params['min_window_size'])
    params['min_particles'] = params['windows'][params['sim']-params['sim_start']][0]
    params['max_particles'] = params['windows'][params['sim']-params['sim_start']][1]

def post_process(params):
    """ Compare the lnpi with the SRSW. """
    num_block = 32
    beta_mu_eq = list()
    rho_vapor = list()
    rho_liquid = list()
    en_vapor = list()
    en_liquid = list()
    pressure = list()
    # convert density of particles/A^3 to kg/m^3 of water
    rho_conv = 1e30/physical_constants.AvogadroConstant().value()*18.01528/1e3
    # convert pressure of kJ/mol/A^3 to bar
    pressure_conv = 1e30/physical_constants.AvogadroConstant().value()*1e3/1e5
    for block in range(-1, num_block):
        if block == -1:
            ln_prob_header = 'ln_prob'
            energy_header = 'e_average'
        else:
            ln_prob_header = 'ln_prob' + str(block)
            energy_header = 'e_block' + str(block)
        lnpi = macrostate_distribution.splice_files(prefix=params['prefix']+'n', suffix='_crit.csv',
                                                    ln_prob_header=ln_prob_header, shift=False)
        lnpi.set_minimum_smoothing(60)
        for node in range(params['num_nodes']):
            prefix = params['prefix']+'n'+str(node)+'s'
            suffix = '_en.csv'
            multistate_accumulator.splice(prefix=prefix, suffix=suffix,
                                          start=node*params['procs_per_node'],
                                          stop=(node+1)*params['procs_per_node'] - 1).to_csv(prefix+suffix)
        energy = multistate_accumulator.splice(prefix=params['prefix']+'n', suffix='s_en.csv',
                                               start=0, stop=params['num_nodes']-1)
        #energy = multistate_accumulator.splice_by_node(prefix=params['prefix']+'n',
        #                                               suffix='_en.txt',
        #                                               num_nodes=params['num_nodes'])
        lnpi.concat_dataframe(dataframe=energy, add_prefix='e_')
        try: # if multiple minima found, then skip the block
            delta_beta_mu = lnpi.equilibrium()
            beta_mu_eq.append(params['beta']*params['mu'] + delta_beta_mu)
            for index, phase in enumerate(lnpi.split()):
                n_gce = phase.average_macrostate()
                rho = rho_conv*n_gce/params['cubic_side_length']**3
                energy = phase.ensemble_average(energy_header)/n_gce
                if index == 0:
                    pressure.append(-pressure_conv*phase.ln_prob()[0]/params['beta']/params['cubic_side_length']**3)
                    rho_vapor.append(rho)
                    en_vapor.append(energy)
                else:
                    rho_liquid.append(rho)
                    en_liquid.append(energy)
        except:
            assert block != -1
    data = pd.DataFrame(data={'rho_vapor': rho_vapor, 'rho_liquid': rho_liquid,
                              'pressure': pressure, 'en_vapor': en_vapor, 'en_liquid': en_liquid,
                              'beta_mu_eq': beta_mu_eq})
    data.to_csv(params['script']+'.csv')
    for col in data.columns:
        print(col, data[col][0], '+/-', data[col][1:].std()/np.sqrt(len(data[col][1:])),
              data[col][1:].mean())

    # skip the following checks if temperature is not 300 K
    if np.abs(params['beta'] - 0.4009078501424201) > 1e-5:
        return

    # check equilibrium properties from https://www.nist.gov/mml/csd/chemical-informatics-group/sat-tmmc-liquid-vapor-coexistence-properties-long-range
    z_factor = 6
    assert np.abs(data['rho_vapor'][0] - 7.373E-03), z_factor*np.sqrt((data['rho_vapor'][1:].std()/np.sqrt(num_block))**2 + 3.253E-05**2)
    assert np.abs(data['rho_liquid'][0] - 9.981E+02), z_factor*np.sqrt((data['rho_liquid'][1:].std()/np.sqrt(num_block))**2 + 2.928E+00**2)
    assert np.abs(data['pressure'][0] - 1.017E-02), z_factor*np.sqrt((data['pressure'][1:].std()/np.sqrt(num_block))**2 + 4.516E-05**2)
    assert np.abs(data['beta_mu_eq'][0] - -1.524E+01), z_factor*np.sqrt((data['beta_mu_eq'][1:].std()/np.sqrt(num_block))**2 + 5.724E-02**2)

    # check lnpi
    z_factor = 15
    lnpi = macrostate_distribution.splice_files(prefix=params['prefix']+'n', suffix='_crit.csv', shift=False)
    srsw = pd.read_csv(params['feasst_install']+'../plugin/flat_histogram/test/data/colMatb0.400908lnz-15.24', skiprows=18, header=None, delim_whitespace=True)
    srsw = pd.concat([lnpi.dataframe(), srsw], axis=1)
    srsw['deltalnPI'] = srsw[1]-srsw[1].shift(1)
    srsw.to_csv(params['prefix']+'_lnpi.csv')
    diverged = srsw[srsw.deltalnPI-srsw.delta_ln_prob > z_factor*srsw.delta_ln_prob_stdev]
    print(diverged)
    assert len(diverged) < 10

    # plot lnpi
    fst = pd.read_csv(params['prefix']+'_lnpi.csv')
    plt.plot(fst['state'], fst['ln_prob'], label='FEASST')
    plt.plot(srsw[0], srsw[1], linestyle='dashed', label='SRSW')
    plt.xlabel('number of particles', fontsize=16)
    plt.ylabel('ln probability', fontsize=16)
    plt.legend(fontsize=16)
    #plt.savefig(params['prefix']+'_lnpi.png', bbox_inches='tight', transparent='True')

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=launch_04_lj_tm_parallel.write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
