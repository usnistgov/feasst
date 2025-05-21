"""
Flat-histogram simulation of patchy trimers in the grand canonical ensemble.
Compare with Fig 9 of http://dx.doi.org/10.1063/1.4918557
This tutorial uses the flat histogram LJ tutorial for most of the simulation setup.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
import launch_04_lj_tm_parallel

def parse():
    """ Parse arguments from command line or change their default values. """
    params, args = launch_04_lj_tm_parallel.parse(
          fstprt='/feasst/particle/trimer.txt',
          beta=1./0.275,
          beta_mu=-5,
          min_sweeps=20,
          cubic_side_length=8,
          max_particles=100,
          window_alpha=1.25,
          procs_per_node=16)
    params['script'] = __file__
    params['prefix'] = 'trimer'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['trial_rigid_cluster_weight'] = 1./params['max_particles']
    params['rwca'] = 2**(1./6.)
    params['system'] = """Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt} cutoff0_1 {rwca} cutoff1_1 {rwca}
WriteModelParams output_file {prefix}_model_params.txt
NeighborCriteria energy_maximum -0.5 site_type0 0 site_type1 0
Potential EnergyMap EnergyMapNeighborCriteria neighbor_index 0 Model LennardJonesForceShift
RefPotential Model LennardJonesForceShift cutoff {rwca} VisitModel VisitModelCell min_length {rwca}""".format(**params)
    params['nvt_trials'] = """TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 0.5 tunable_param 0.2 tunable_target_acceptance 0.25 particle_type 0
TrialRigidCluster weight {trial_rigid_cluster_weight} neighbor_index 0""".format(**params)
    params['muvt_trials'] = "TrialTransfer weight 2 particle_type 0 reference_index 0 num_steps 8"
    return params, args

def post_process(params):
    """ Test not implemented """
    assert True

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=launch_04_lj_tm_parallel.sim_node_dependent_params,
                          write_feasst_script=launch_04_lj_tm_parallel.write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
