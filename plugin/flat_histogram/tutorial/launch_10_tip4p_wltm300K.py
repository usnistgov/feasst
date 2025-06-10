"""
This tutorial is similar to tutorial 10 spce, but for tip4p.
"""

import argparse
from pyfeasst import fstio
from pyfeasst import physical_constants
from pyfeasst import macrostate_distribution
import launch_04_lj_tm_parallel
import launch_10_spce_wltm300K

def parse(temperature=300):
    """ Parse arguments from command line or change their default values. """
    params, args = launch_04_lj_tm_parallel.parse(
          fstprt='/feasst/plugin/charge/particle/tip4p.txt',
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
    params['prefix'] = 'tip4p_lowt'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['dccb_cut'] = 0.75*3.165
    params['ewald_alpha'] = 5.6/params['cubic_side_length']
    params['num_particles_first_node'] = 180

    # write TrialGrowFile
    with open(params['prefix']+'_grow.txt', 'w') as f:
        f.write("""TrialGrowFile

particle_type=water weight=2 transfer=true site=O1 num_steps=10 reference_index=0
bond=true mobile_site=M1 anchor_site=O1 reference_index=0
branch=true mobile_site2=H1 mobile_site=H2 anchor_site=O1 anchor_site2=M1 reference_index=0
""")

    return params, args

def post_process(params):
    import numpy as np
    lnpi = macrostate_distribution.splice_files(prefix=params['prefix']+'n', suffix='_crit.csv', shift=False)
    #lnpi.plot(show=True)
    lnpi.set_minimum_smoothing(50)
    lnpi.reweight(delta_beta_mu=1.5, inplace=True)
    lnpi.equilibrium()
    #lnpi.plot(show=True)
    vapor, liquid = lnpi.split()
    print('<N_vapor>_GCE', vapor.average_macrostate())
    assert np.abs(vapor.average_macrostate() - 0.01) < 0.001
    print('<N_liquid>_GCE', liquid.average_macrostate())
    assert np.abs(liquid.average_macrostate() - 264.1) < 1.5

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=launch_10_spce_wltm300K.sim_node_dependent_params,
                          write_feasst_script=launch_04_lj_tm_parallel.write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
