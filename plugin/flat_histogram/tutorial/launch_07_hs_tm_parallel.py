"""
Flat-histogram simulation of single-site hard sphere particles in the grand canonical ensemble.
Compare equation of state in SRSW:
https://www.nist.gov/mml/csd/chemical-informatics-research-group/hard-sphere-thermodynamic-and-transport-properties
This tutorial uses flat histogram LJ tutorial for most of the simulation setup.
"""

import argparse
import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import macrostate_distribution
import launch_04_lj_tm_parallel

def parse():
    params, args = launch_04_lj_tm_parallel.parse(
          fstprt='/feasst/particle/atom.txt',
          beta=1.,
          beta_mu=-2.352321,
          min_sweeps=1e2,
          cubic_side_length=8,
          max_particles=256,
          window_alpha=1.5)
    """ Parse arguments from command line or change their default values. """
    params['script'] = __file__
    params['prefix'] = 'hs'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['system'] = """Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt}
Potential Model HardSphere VisitModel VisitModelCell min_length 1""".format(**params)
    params['nvt_trials'] = "TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25"
    return params, args

def post_process(params):
    # compare to EOS in SRSW: https://www.nist.gov/mml/csd/chemical-informatics-research-group/hard-sphere-thermodynamic-and-transport-properties
    lnpi=macrostate_distribution.splice_files(prefix=params['prefix']+'n0s', suffix='_crit.csv', shift=False)
    #lnpi = macrostate_distribution.MacrostateDistribution(file_name=params['prefix']+'n0_lnpi.txt')
    volume = 8**3
    srsw = pd.read_csv(params['feasst_install']+'../plugin/flat_histogram/test/data/stat_hs.csv')
    srsw = srsw[:6]
    pressure = list()
    lnpi_rw = copy.deepcopy(lnpi)
    for target_density in srsw['dens']:
        lnpi_rw.reweight_to_macrostate(target_macrostate=target_density*volume)
        pressure.append(-lnpi_rw.ln_prob().values[0]/volume)
    srsw['P_FST'] = pressure
    print(srsw[['dens', 'P_MC', 'P_FST', '+/-']])
    assert np.any(abs(srsw['P_MC'] - srsw['P_FST']) < 1e-3)

    # Use chemical potential from Carnahan-Starling to compare expected average density
    # http://www.sklogwiki.org/SklogWiki/index.php/Carnahan-Starling_equation_of_state
    rho = 0.1
    cubic_side_length = 8
    eta = np.pi/6*rho
    betamu_ex = (8*eta-9*eta**2+3*eta**3)/(1-eta)**3
    betamu = betamu_ex + np.log(rho)
    lnpi_rw = lnpi.reweight(delta_beta_mu=betamu+2.352321)
    density = lnpi_rw.average_macrostate()/volume
    print('target_density', rho, 'density', density)
    assert abs(rho - density) < 5e-4

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=launch_04_lj_tm_parallel.sim_node_dependent_params,
                          write_feasst_script=launch_04_lj_tm_parallel.write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
