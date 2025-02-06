"""
Flat-histogram simulation of EMP2 CO2 (https://doi.org/10.1021/j100031a034) in the grand canonical ensemble.
This tutorial uses the flat histogram LJ tutorial for most of the simulation setup.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import physical_constants
from pyfeasst import macrostate_distribution
import launch_04_lj_tm_parallel

def parse(temperature=298):
    """ Parse arguments from command line or change their default values. """
    params, args = launch_04_lj_tm_parallel.parse(
          fstprt='/feasst/particle/co2_epm2.fstprt',
          beta=1./(temperature*physical_constants.MolarGasConstant().value()/1e3), # mol/kJ
          beta_mu=-6,
          min_sweeps=2,
          cubic_side_length=28,
          max_particles=300,
          window_alpha=2.15,
          collect_flatness=18,
          min_flatness=22,
          hours_checkpoint=0.2,
          hours_terminate=1,
          min_window_size=3)
    params['script'] = __file__
    params['prefix'] = 'co2'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['min_particles_second_window'] = "min1 40"
    params['cutoff'] = 12
    params['dccb_cut'] = 4.
    params['dccb_cut'] = params['cubic_side_length']/int(params['cubic_side_length']/params['dccb_cut']) # maximize inside box
    params['ewald_alpha'] = 5.6/params['cubic_side_length']
    params['mu_init']=10
    params['angle_width'] = 0.01
    params['angle_center'] = np.pi - 0.5*params['angle_width']
    params['system'] = """Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt} cutoff {cutoff} sigma0_1 2.89170901025674 group0 oxygen oxygen_site_type 0
Potential VisitModel Ewald alpha {ewald_alpha} kmax_squared 38
Potential Model ModelTwoBodyFactory model0 LennardJones model1 ChargeScreened erfc_table_size 2e4 VisitModel VisitModelCutoffOuter
RefPotential VisitModel DontVisitModel
#RefPotential Model HardSphere group oxygen cutoff {dccb_cut} VisitModel VisitModelCell min_length {dccb_cut} cell_group oxygen
Potential Model ChargeScreenedIntra VisitModel VisitModelBond
Potential Model ChargeSelf
Potential VisitModel LongRangeCorrections""".format(**params)
    params['nvt_trials'] = """TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 0.25 particle_type 0 tunable_param 0.2 tunable_target_acceptance 0.25
TrialGrowFile grow_file {prefix}_grow_canonical.txt""".format(**params)
    params['init_trials'] = """TrialGrowFile grow_file {prefix}_grow_grand_canonical.txt""".format(**params)
    params['init_remove'] = "Remove name_contains0 add name_contains1 remove"
    params['muvt_trials'] = """TrialGrowFile grow_file {prefix}_grow_grand_canonical.txt
#AnalyzeBonds trials_per_write {tpc} output_file {prefix}n{node}s[sim_index]_bonds.txt start_after_cycle 1 angle_bin_center {angle_center} angle_bin_width {angle_width} bond_bin_center 1.149 bond_bin_width 0.00001
#PairDistribution trials_per_update 1000 trials_per_write {tpc} start_after_cycle 1 dr 0.025 output_file {prefix}n{node}s[sim_index]_gr.txt multistate true multistate_aggregate false""".format(**params)

    # write TrialGrowFile
    with open(params['prefix']+'_grow_grand_canonical.txt', 'w') as f:
        f.write("""TrialGrowFile

particle_type 0 weight 2 transfer true site 0 num_steps 1 reference_index 0
bond true mobile_site 1 anchor_site 0 num_steps 1 reference_index 0
angle true mobile_site 2 anchor_site 0 anchor_site2 1 reference_index 0
""")
    with open(params['prefix']+'_grow_canonical.txt', 'w') as f:
        f.write("""TrialGrowFile

particle_type 0 weight 0.1 angle true mobile_site 2 anchor_site 0 anchor_site2 1 reference_index 0

particle_type 0 weight 0.1 angle true mobile_site 1 anchor_site 0 anchor_site2 2 reference_index 0
""")

    return params, args

def post_process(params):
    lnp = macrostate_distribution.splice_files(prefix=params['prefix']+'n0s', suffix='_crit.csv', shift=False)
#    lnp.reweight(-0.75, inplace=True)
#    delta_beta_mu = lnp.equilibrium(delta_beta_mu_guess=0.01)
#    #print(delta_beta_mu)
#    #print(lnp.dataframe())
#    #lnp.plot()
#    #plt.savefig(params['prefix']+'.png')
#    #print(lnp.minimums())
#    vapor, liquid = lnp.split()
#    volume = params['cubic_side_length']**3
#    na = physical_constants.AvogadroConstant().value()
#    dens_conv = 1./volume/na*44.01/1e3*1e30 # convert from N/V units of molecules/A^3 to kg/m^3
#    print('vapor density(kg/m^3)', vapor.average_macrostate()*dens_conv)
#    print('liquid density(kg/m^3)', liquid.average_macrostate()*dens_conv)
#    assert np.abs(256 - vapor.average_macrostate()*dens_conv) < 10
#    assert np.abs(712 - liquid.average_macrostate()*dens_conv) < 50

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=launch_04_lj_tm_parallel.sim_node_dependent_params,
                          write_feasst_script=launch_04_lj_tm_parallel.write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
