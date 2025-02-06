"""
This tutorial is similar to tutorial 4, but for this low temperature simulation,
we will split the simulation into two different nodes.
The first node will have less particles but a higher number of sweeps required.
The second node will have dccb but not avb.
This tutorial uses flat histogram LJ tutorial for most of the simulation setup.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import macrostate_distribution
from pyfeasst import multistate_accumulator
import launch_04_lj_tm_parallel

def parse():
    """ Parse arguments from command line or change their default values. """
    params, args = launch_04_lj_tm_parallel.parse(
          fstprt='/feasst/particle/lj.fstprt',
          beta=1./0.7,
          beta_mu=-5.943376,
          min_sweeps=5,
          cubic_side_length=8,
          max_particles=475,
          window_alpha=1.5,
          hours_checkpoint=0.1,
          hours_terminate=5*24,
          num_nodes=2,
          collect_flatness=18,
          min_flatness=22)
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    params['script'] = __file__
    params['prefix'] = 'lj_lowt'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['dccb_cut_min'] = 1
    params['dccb_cut'] = params['cubic_side_length']/int(params['cubic_side_length']/params['dccb_cut_min'])
    params['num_particles_first_node'] = 375
    return params, args

def sim_node_dependent_params(params):
    """ Define parameters that are dependent on the sim or node. """
    if params['node'] == 0:
        params['min_particles'] = 0
        params['max_particles'] = params['num_particles_first_node']
        params['muvt_trials'] = 'TrialTransfer weight 2 particle_type 0\nTrialTransferAVB weight 0.2 particle_type 0'
        params['lj_potential'] = 'Potential EnergyMap EnergyMapNeighborCriteria neighbor_index 0 Model LennardJones'
        params['ref_potential'] = ''
        params['avb_trials'] = 'TrialAVB2 weight 0.1 particle_type 0\nTrialAVB4 weight 0.1 particle_type 0'
        params['min_sweeps'] = 20
        params['window_alpha'] = 2
        params['min_window_size'] = 5
    elif params['node'] == 1:
        params['min_particles'] = params['num_particles_first_node']
        params['muvt_trials'] = 'TrialTransfer weight 2 particle_type 0 reference_index 0 num_steps 10'
        params['lj_potential'] = 'Potential Model LennardJones'
        params['ref_potential'] = """RefPotential Model LennardJones cutoff {dccb_cut} VisitModel VisitModelCell min_length {dccb_cut}""".format(**params)
        params['avb_trials'] = ''
        params['min_sweeps'] = 2
        params['window_alpha'] = 1
        params['min_window_size'] = 3
    params['system'] = """Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt}
NeighborCriteria maximum_distance 1.375 minimum_distance 0.9
{lj_potential}
{ref_potential}
Potential VisitModel LongRangeCorrections""".format(**params)
    params['nvt_trials'] = """TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
{avb_trials}""".format(**params)
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
    for block in range(-1, num_block):
        if block == -1:
            ln_prob_header = 'ln_prob'
            energy_header = 'e_average'
        else:
            ln_prob_header = 'ln_prob' + str(block)
            energy_header = 'e_block' + str(block)
        lnpi = macrostate_distribution.splice_files(prefix=params['prefix']+'n', suffix='_crit.csv',
                                                    ln_prob_header=ln_prob_header, shift=False)
        lnpi.set_minimum_smoothing(30)
        for node in range(params['num_nodes']):
            prefix = params['prefix']+'n'+str(node)+'s'
            suffix = '_en.csv'
            multistate_accumulator.splice(prefix=prefix, suffix=suffix,
                                          start=node*params['procs_per_node'],
                                          stop=(node+1)*params['procs_per_node'] - 1).to_csv(prefix+suffix)
        energy = multistate_accumulator.splice(prefix=params['prefix']+'n', suffix='s_en.csv',
                                               start=0, stop=params['num_nodes']-1)
        lnpi.concat_dataframe(dataframe=energy, add_prefix='e_')
        delta_beta_mu = lnpi.equilibrium()
        beta_mu_eq.append(params['beta']*params['mu'] + delta_beta_mu)
        for index, phase in enumerate(lnpi.split()):
            n_gce = phase.average_macrostate()
            rho = n_gce/params['cubic_side_length']**3
            energy = phase.ensemble_average(energy_header)/n_gce
            if index == 0:
                pressure.append(-phase.ln_prob()[0]/params['beta']/params['cubic_side_length']**3)
                rho_vapor.append(rho)
                en_vapor.append(energy)
            else:
                rho_liquid.append(rho)
                en_liquid.append(energy)
    data = pd.DataFrame(data={'rho_vapor': rho_vapor, 'rho_liquid': rho_liquid,
                              'pressure': pressure, 'en_vapor': en_vapor, 'en_liquid': en_liquid,
                              'beta_mu_eq': beta_mu_eq})
    data.to_csv(params['script']+'.csv')
    for col in data.columns:
        print(col, data[col][0], '+/-', data[col][1:].std()/np.sqrt(len(data[col][1:])),
              data[col][1:].mean())

    # skip the following checks if temperature is not 0.7
    if np.abs(params['beta'] - 1./0.7) > 1e-5:
        return

    # check equilibrium properties from https://www.nist.gov/mml/csd/chemical-informatics-group/sat-tmmc-liquid-vapor-coexistence-properties-long-range
    z_factor = 6
    assert np.abs(data['rho_vapor'][0] - 1.996E-03) < z_factor*np.sqrt((data['rho_vapor'][1:].std()/np.sqrt(num_block))**2 + 1.422E-05**2)
    assert np.abs(data['rho_liquid'][0] - 8.437E-01), z_factor*np.sqrt((data['rho_liquid'][1:].std()/np.sqrt(num_block))**2 + 2.49E-04**2)
    assert np.abs(data['pressure'][0] - 1.370E-03), z_factor*np.sqrt((data['pressure'][1:].std()/np.sqrt(num_block))**2 + 9.507E-07**2)
    assert np.abs(data['en_vapor'][0] - -2.500E-02), z_factor*np.sqrt((data['en_vapor'][1:].std()/np.sqrt(num_block))**2 + 3.323E-05**2)
    assert np.abs(data['en_liquid'][0] - -6.106E+00), z_factor*np.sqrt((data['en_liquid'][1:].std()/np.sqrt(num_block))**2 + 1.832E-03**2)
    assert np.abs(data['beta_mu_eq'][0] - -6.257E+00), z_factor*np.sqrt((data['beta_mu_eq'][1:].std()/np.sqrt(num_block))**2 + 4.639E-04**2)

    # check lnpi
    z_factor = 25
    lnpi = macrostate_distribution.splice_files(prefix=params['prefix']+'n', suffix='_crit.csv', shift=False)
    srsw = pd.read_csv(params['feasst_install']+'../plugin/flat_histogram/test/data/stat070.csv')
    srsw = pd.concat([lnpi.dataframe(), srsw], axis=1)
    srsw['deltalnPI'] = srsw.lnPI-srsw.lnPI.shift(1)
    srsw.to_csv(params['prefix']+'_lnpi.csv')
    diverged = srsw[srsw.deltalnPI-srsw.delta_ln_prob > z_factor*srsw.delta_ln_prob_stdev]
    print('diverged', diverged)
    assert len(diverged) < 10

    # plot lnpi
    fst = pd.read_csv(params['prefix']+'_lnpi.csv')
    plt.plot(fst['state'], fst['ln_prob'], label='FEASST')
    plt.plot(srsw['N'], srsw['lnPI'], linestyle='dashed', label='SRSW')
    plt.xlabel('number of particles', fontsize=16)
    plt.ylabel('ln probability', fontsize=16)
    plt.legend(fontsize=16)
    plt.savefig(params['prefix']+'_lnpi.png', bbox_inches='tight', transparent='True')

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=launch_04_lj_tm_parallel.write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
