"""
Flat-histogram simulation of TraPPE alkanes in the grand canonical ensemble.
https://www.nist.gov/mml/csd/chemical-informatics-group/sat-tmmc-liquid-vapor-coexistence-properties-trappe-ua-n-butane
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

def parse(temperature=350):
    """ Parse arguments from command line or change their default values. """
    params, args = launch_04_lj_tm_parallel.parse(
          fstprt='/feasst/particle/n-octane.fstprt',
          beta=1./(temperature*physical_constants.MolarGasConstant().value()/1e3), # mol/kJ
          beta_mu=-6,
          min_sweeps=2,
          cubic_side_length=45,
          max_particles=285,
          window_alpha=1.5,
          collect_flatness=18,
          min_flatness=22,
          hours_checkpoint=1,
          hours_terminate=5*24,
          min_window_size=3)
    params['script'] = __file__
    params['prefix'] = 'trappe'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['windows'] = macrostate_distribution.window_exponential(
        alpha=params['window_alpha'], minimums=[params['min_particles'], 17], maximum=params['max_particles'],
        number=params['num_sims'], overlap=1, min_size=params['min_window_size'])
    params['cutoff'] = 12
    params['dccb_cut'] = 4.
    params['dccb_cut'] = params['cubic_side_length']/int(params['cubic_side_length']/params['dccb_cut']) # maximize inside box
    params['num_sites'] = fstio.num_sites_in_fstprt(params['fstprt'], params['feasst_install'])
    if 'n-butane' in params['fstprt']:
        params['molecular_weight'] = 58.12
    elif 'n-octane' in params['fstprt']:
        params['molecular_weight'] = 114.23
    else:
        assert False, "input new molecular_weight into params"
    params['last_site'] = params['num_sites'] - 1
    params['system'] = """Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt} cutoff {cutoff}
Potential Model LennardJones
Potential Model LennardJones VisitModel VisitModelIntra intra_cut 3
Potential VisitModel LongRangeCorrections
RefPotential Model LennardJones VisitModel VisitModelCell min_length {dccb_cut} reference_index 0
RefPotential Model LennardJones VisitModel VisitModelIntra intra_cut 3 reference_index 0""".format(**params)
    params['nvt_trials'] = """TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 0.25 particle_type 0 tunable_param 0.2 tunable_target_acceptance 0.25 pivot_site 0
TrialParticlePivot weight 0.25 particle_type 0 tunable_param 0.2 tunable_target_acceptance 0.25 pivot_site {last_site}
TrialGrowFile grow_file trappe_grow_canonical.txt""".format(**params)
    params['init_trials'] = "TrialGrowFile grow_file trappe_grow_grand_canonical.txt"
    params['init_remove'] = "Remove name_contains0 add name_contains1 remove"
    params['muvt_trials'] = "TrialGrowFile grow_file trappe_grow_grand_canonical.txt"
    write_grow_file(filename="trappe_grow_canonical.txt", params=params, gce=False)
    write_grow_file(filename="trappe_grow_grand_canonical.txt", params=params, gce=True)
    return params, args

def write_partial(f, bond, angle, dihedral, params):
    if params['num_sites'] == 2:
        f.write(bond)
    elif params['num_sites'] == 3:
        f.write(angle)
    elif params['num_sites'] > 3:
        f.write(dihedral)
    else:
        print('unrecognized num_sites', params['num_sites'])
        assert False

# write TrialGrowFile to include grand canonical ensemble growth and canonica ensemble reptations
def write_grow_file(filename, params, gce):
    with open(filename, 'w') as f:
        f.write("TrialGrowFile\n\n")
        for inv in [True, False]:
            for trial_type in range(3+int(params['num_sites']/2)): # 0: reptate, 1: full regrow, 2+: partial regrow
            #for trial_type in [0, 1, 2]: # 0: reptate, 1: full regrow, 2: partial regrow
                for site in range(params['num_sites']):
                    for i in range(4):
                        sign = -1
                        #if (trial_type == 0 or trial_type == 2) and site != params['num_sites'] - 1:
                        if trial_type == 0 and site != params['num_sites'] - 1:
                            sign = 1
                        params['site'+str(i)] = site + sign*i
                        if inv:
                            params['site'+str(i)] = params['num_sites'] - site - 1 - sign*i
                    bond = """bond true mobile_site {site0} anchor_site {site1} num_steps 4 reference_index 0\n""".format(**params)
                    angle = """angle true mobile_site {site0} anchor_site {site1} anchor_site2 {site2} num_steps 4 reference_index 0\n""".format(**params)
                    dihedral = """dihedral true mobile_site {site0} anchor_site {site1} anchor_site2 {site2} anchor_site3 {site3} num_steps 4 reference_index 0\n""".format(**params)

                    # full regrowth insertion/deletion
                    if trial_type == 1 and gce:
                        if site == 0:
                            f.write("""particle_type 0 weight 2 transfer true site {site0} num_steps 4 reference_index 0\n""".format(**params))
                        elif site == 1:
                            f.write(bond)
                        elif site == 2:
                            f.write(angle)
                        else:
                            f.write(dihedral)

#                    # reptation. There seems to be a problem with reptation.
#                    elif trial_type == 0 and not gce:
#                        if site == params['num_sites'] - 1:
#                            write_partial(f, bond, angle, dihedral, params)
#                        else:
#                            if site == 0:
#                                f.write("""particle_type 0 weight 2 """)
#                            f.write("""reptate true mobile_site {site0} anchor_site {site1} num_steps 1 reference_index 0\n""".format(**params))

                    # partial regrow
                    if not gce and trial_type > 1:
                        num_grow = trial_type - 1
                        if params['num_sites'] - site < num_grow:
                            if params['num_sites'] - site == num_grow - 1:
                                f.write('particle_type 0 weight '+str(2/(trial_type-2))+' ')
                            if site == 1:
                                f.write(bond)
                            elif site == 2:
                                f.write(angle)
                            elif site != 0:
                                f.write(dihedral)

                f.write("\n")

def post_process(params):
    #lnp = macrostate_distribution.splice_collection_matrix(prefix='trappen0s', suffix='_crit.txt', use_soft=True)
    lnp = macrostate_distribution.splice_files(prefix=params['prefix']+'n0s', suffix='_crit.csv', shift=False)
    assert np.abs(lnp.ln_prob()[1] - lnp.ln_prob()[0] - 5.86440399999992) < 0.1
    assert np.abs(lnp.ln_prob()[2] - lnp.ln_prob()[1] - 5.22683300000017) < 0.1
    # equilibrium test below was abandoned to reduce max_particles for faster convergence
    #lnp.equilibrium()
    #lnp.plot(show=True)
    #print('WARNING: max_particles should be higher but the liquid peak was truncated to make the simulation faster')
    #vapor, liquid = lnp.split()
    #volume = params['cubic_side_length']**3
    #na = physical_constants.AvogadroConstant().value()
    #dens_conv = 1./volume/na*params['molecular_weight']/1e3*1e30 # convert from N/V units of molecules/A^3 to kg/m
    ## https://www.nist.gov/mml/csd/chemical-informatics-research-group/sat-tmmc-liquid-vapor-coexistence-properties-trappe-ua-n
    #density = vapor.average_macrostate()*dens_conv
    #print('density', density)
    #assert np.abs(30.6 - density) < 2
    #assert np.abs(508 - liquid.average_macrostate()*dens_conv) < 30

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=launch_04_lj_tm_parallel.sim_node_dependent_params,
                          write_feasst_script=launch_04_lj_tm_parallel.write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
