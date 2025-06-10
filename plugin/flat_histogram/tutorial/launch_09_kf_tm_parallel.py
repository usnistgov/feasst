"""
Flat-histogram simulation of Kern-Frenkel patchy particles in the grand canonical ensemble.
https://doi.org/10.1063/1.1569473
This tutorial uses the flat histogram LJ tutorial for most of the simulation setup.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import macrostate_distribution
import launch_04_lj_tm_parallel

def parse():
    """ Parse arguments from command line or change their default values. """
    params, args = launch_04_lj_tm_parallel.parse(
          fstprt='/feasst/plugin/patch/particle/two_patch_linear.txt',
          beta=1./0.7,
          beta_mu=-2.14285714285714,
          min_sweeps=5,
          cubic_side_length=8,
          max_particles=370,
          window_alpha=2.25)
    params['script'] = __file__
    params['prefix'] = 'kf'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['chi'] = 0.7
    params['patch_angle'] = 2*np.arcsin(np.sqrt(params['chi']/2))*180/np.pi
    params['system'] = """Configuration cubic_side_length={cubic_side_length} particle_type=lj:{fstprt} \
  patch_angleP={patch_angle} group=centers centers_site_type=C
Potential Model=HardSphere VisitModel=VisitModelCell min_length=1 cell_group=centers group=centers
Potential Model=SquareWell VisitModel=VisitModelCell min_length=1.5 cell_group=centers \
  VisitModelInner=VisitModelInnerPatch group=centers""".format(**params)
    params['nvt_trials'] = """TrialTranslate weight=1 tunable_param=0.2
TrialRotate weight=1 tunable_param=0.2"""
    return params, args

def post_process(params):
    lnpi=macrostate_distribution.splice_files(prefix=params['prefix']+'n0s', suffix='_crit.csv', shift=False)
    lnpi=lnpi.dataframe()
    lnpi = lnpi[:6] # cut down to six rows
    lnpi['ln_prob'] -= np.log(sum(np.exp(lnpi['ln_prob'])))  # renormalize
    lnpi['ln_prob_prev'] = [-15.9976474469475, -11.9104563420586,  -8.48324267323538, -5.42988602574393, -2.64984051640555, -0.07824246342703]
    lnpi['ln_prob_prev_stdev'] = [0.035, 0.03, 0.025, 0.02, 0.015, 0.01]
    diverged = lnpi[lnpi.ln_prob-lnpi.ln_prob_prev > 6*lnpi.ln_prob_prev_stdev]
    print(diverged)
    assert len(diverged) == 0
    energy = pd.read_csv(params['prefix']+'n0s000_en.csv')
    energy = energy[:6]
    energy['prev'] = [0, 0, -0.038758392176564, -0.116517384264731, -0.232665619265520, -0.387804181572135]
    diverged = energy[energy.average - energy.prev > 10*energy.block_stdev]
    print(diverged)
    assert len(diverged) == 0

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=launch_04_lj_tm_parallel.sim_node_dependent_params,
                          write_feasst_script=launch_04_lj_tm_parallel.write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
