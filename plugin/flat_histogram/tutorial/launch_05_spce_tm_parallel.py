"""
Flat-histogram simulation of SPC/E water in the grand canonical ensemble.
The default temperature of 525 K is above the critical point.
The results are are compared with the SRSW https://doi.org/10.18434/T4M88Q
(which use CODATA2010 physical constants).
This tutorial uses flat histogram LJ tutorial for most of the simulation setup.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import physical_constants
from pyfeasst import macrostate_distribution
import launch_04_lj_tm_parallel

def parse(temperature=525):
    """ Parse arguments from command line or change their default values. """
    params, args = launch_04_lj_tm_parallel.parse(
          fstprt='/feasst/particle/spce.txt',
          beta=1./(temperature*physical_constants.MolarGasConstant().value()/1e3), # mol/kJ
          beta_mu=-8.14,
          min_sweeps=5,
          cubic_side_length=20,
          max_particles=265,
          window_alpha=1.5)
    params['script'] = __file__
    params['prefix'] = 'spce'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['ewald_alpha'] = 5.6/params['cubic_side_length']
    params['system'] = """Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt} group0 oxygen oxygen_site_type 0
Potential VisitModel Ewald alpha {ewald_alpha} kmax_squared 38
Potential Model ModelTwoBodyFactory model0 LennardJones model1 ChargeScreened erfc_table_size 2e4 VisitModel VisitModelCutoffOuter
Potential Model ChargeScreenedIntra VisitModel VisitModelBond
Potential Model ChargeSelf
Potential VisitModel LongRangeCorrections""".format(**params)
    params['nvt_trials'] = """TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 0.5 particle_type 0 tunable_param 0.5 tunable_target_acceptance 0.25"""
    params['muvt_trials'] = """TrialTransfer weight 2 particle_type 0"""
    return params, args

def post_process(params):
    """ Skip the following checks if temperature is not 525 K """
    if np.abs(params['beta'] - 0.22909020008138298) > 1e-5:
        return
    #lnpi = pd.read_csv(params['prefix']+'n0_lnpi.txt')
    lnpi=macrostate_distribution.splice_files(prefix=params['prefix']+'n0s', suffix='_crit.csv', shift=False)
    lnpi=lnpi.dataframe()
    lnpi = pd.concat([lnpi, pd.read_csv(params['feasst_install']+'../plugin/flat_histogram/test/data/stat_spce_525.csv')], axis=1)
    lnpi['deltalnPI'] = lnpi.lnPI - lnpi.lnPI.shift(1)
    diverged = lnpi[lnpi.deltalnPI-lnpi.delta_ln_prob > 6*lnpi.delta_ln_prob_stdev]
    print(diverged)
    assert len(diverged) < int(0.15*len(lnpi))
    plt.plot(lnpi['state'], lnpi['ln_prob'], label='FEASST')
    plt.plot(lnpi['N'], lnpi['lnPI'], linestyle='dashed', label='SRSW')
    plt.xlabel('number of particles', fontsize=16)
    plt.ylabel('ln probability', fontsize=16)
    plt.legend(fontsize=16)
    #plt.savefig(params['prefix']+'_lnpi.png', bbox_inches='tight', transparent='True')

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=launch_04_lj_tm_parallel.sim_node_dependent_params,
                          write_feasst_script=launch_04_lj_tm_parallel.write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
