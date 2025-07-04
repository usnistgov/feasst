"""
Flat-histogram simulation of a binary mixture of CO2 and N2 using SAFT-based MIE models in the grand canonical ensemble.
https://doi.org/10.1021/acs.jpcb.5c00536
https://doi.org/10.1063/1.1844372
https://doi.org/10.1063/1.2064628
https://doi.org/10.1063/1.4966573
https://doi.org/10.1063/1.4996759
https://doi.org/10.1063/1.5006906

In the postprocessing step, reweight in mu1 to find equilibrium between vapor and liquid.
The chosen delta mu = mu2 - mu1 determines the resulting pressure.
Simulations of various delta mu yeild the binary phase envelope for a given temperature.

Each macrostate is an isochoric semigrand ensemble (instead of canonical in the single component case), but we are simulating with grand canonical ensemble sampling (e.g., insertions and deletions).
Note that when computing mole fractions, <N_0>/<N_total> != <N_0/N_total>, where <...> is a isochoric semigrand ensemble average.
Use gce ensemble averages of extensive quantities, not intensive quantities.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import physical_constants
from pyfeasst import macrostate_distribution
from pyfeasst import multistate_accumulator

def parse(particle_type1='/feasst/particle/dimer_mie_CO2.txt',
          particle_type2='/feasst/particle/dimer_mie_N2.txt',
          beta=1./258.15,
          beta_mu1=-3,
          beta_delta_mu=0.244046,
          beta_init=1/400,
          min_sweeps=5,
          cubic_side_length=24,
          max_particles=215,
          window_alpha=2.25,
          hours_checkpoint=1,
          hours_terminate=5*24,
          num_nodes=1,
          min_window_size=3,
          procs_per_node=32,
          tpc=int(1e6),
          collect_flatness=20,
          min_flatness=25):
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--particle_type1', type=str, default=particle_type1, help='first component FEASST particle')
    parser.add_argument('--particle_type2', type=str, default=particle_type2, help='second component FEASST particle')
    parser.add_argument('--beta', type=float, default=beta, help='1/(Temperature in Kelvin)')
    parser.add_argument('--beta_mu1', type=float, default=beta_mu1, help='chemical potential of first component')
    parser.add_argument('--beta_delta_mu', type=float, default=beta_delta_mu, help='chemical potential of second component minus first')
    parser.add_argument('--beta_init', type=float, default=beta_init, help='chemical potential to initialize number of particles')
    parser.add_argument('--max_particles', type=int, default=max_particles, help='maximum number of particles')
    parser.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
    parser.add_argument('--window_alpha', type=int, default=window_alpha, help='as window_alpha increases, window size of large N decreases')
    parser.add_argument('--min_window_size', type=int, default=min_window_size, help='minimum window size when parallelizing')
    parser.add_argument('--min_sweeps', type=int, default=min_sweeps,
                        help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
    parser.add_argument('--collect_flatness', type=int, default=collect_flatness, help='WL flatness to begin collection matrix')
    parser.add_argument('--min_flatness', type=int, default=min_flatness, help='WL flatness to begin transition matrix')
    parser.add_argument('--cubic_side_length', type=float, default=cubic_side_length,
                        help='cubic periodic boundary length')
    parser.add_argument('--tpc', type=int, default=tpc, help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=0, help='number of cycles for equilibration')
    parser.add_argument('--hours_checkpoint', type=float, default=hours_checkpoint, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=hours_terminate, help='hours until termination')
    parser.add_argument('--procs_per_node', type=int, default=procs_per_node, help='number of processors')
    parser.add_argument('--run_type', '-r', type=int, default=0,
                        help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--seed', type=int, default=-1,
                        help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
    parser.add_argument('--num_nodes', type=int, default=num_nodes, help='Number of nodes in queue')
    parser.add_argument('--scratch', type=str, default=None,
                        help='Optionally write scheduled job to scratch/logname/jobid.')
    parser.add_argument('--node', type=int, default=0, help='node ID')
    parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
    parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

    # Convert arguments into a parameter dictionary, and add argument-dependent parameters.
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    params = vars(args)
    params['script'] = __file__
    params['prefix'] = 'binary'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.95*params['hours_terminate'] - 0.05 # terminate FEASST before SLURM
    params['procs_per_sim'] = 1
    params['num_sims'] = params['procs_per_node']*params['num_nodes']
    params['mu1'] = params['beta_mu1']/params['beta']
    params['mu2'] = params['beta_delta_mu']/params['beta'] + params['mu1']
    params['min_particles_second_window'] = ""
    params['minimums'] = [params['min_particles']]
    if params['num_nodes'] == 1:
        params['windows'] = macrostate_distribution.window_exponential(
            alpha=params['window_alpha'], minimums=params['minimums'], maximum=params['max_particles'],
            number=params['num_sims'], overlap=1, min_size=params['min_window_size'])
    return params, args

def sim_node_dependent_params(params):
    """ Define parameters that are dependent on the sim or node. """
    params['min_particles'] = params['windows'][params['sim']][0]
    params['max_particles'] = params['windows'][params['sim']][1]
    params['sim_start'] = 0
    params['sim_end'] = params['num_sims'] - 1

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed={seed}
Configuration cubic_side_length={cubic_side_length} particle_type=pt1:{particle_type1},pt2:{particle_type2} \
    model_param_file=/feasst/particle/mie_model_parameters.txt
WriteModelParams output_file={prefix}{sim:03d}_model_params.txt
Potential Model=Mie table_size=1e4
Potential VisitModel=LongRangeCorrections
RefPotential VisitModel=DontVisitModel ref=noixn
ThermoParams beta={beta_init} chemical_potential={mu1},{mu2}
Metropolis
For [pt]=pt1,pt2
    TrialTranslate weight_per_number_fraction=0.5 particle_type=[pt] tunable_param=0.2
    TrialParticlePivot weight_per_number_fraction=0.5 particle_type=[pt] tunable_param=0.5
EndFor
# The following semigrand particle type changes only help sampling if the excluded volumes of the two particles are similar (and particles are rigid)
TrialMorph weight=0.5 particle_type=pt1 particle_type_morph=pt2 ref=noixn
TrialMorph weight=0.5 particle_type=pt2 particle_type_morph=pt1 ref=noixn
CheckEnergy trials_per_update={tpc} decimal_places=4
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# gcmc initialization and nvt equilibration
TrialAdd particle_type=pt1
TrialAdd particle_type=pt2
Let [write]=trials_per_write={tpc} output_file={prefix}n{node}s{sim:03d}
Log [write]_eq.csv
Tune
Run until_num_particles={min_particles}
Remove name=TrialAdd,TrialAdd
ThermoParams beta={beta} chemical_potential={mu1},{mu2}
Metropolis trials_per_cycle={tpc} cycles_to_complete={equilibration}
Run until=complete
Remove name=Tune,Log

# gcmc tm production
FlatHistogram Macrostate=MacrostateNumParticles width=1 max={max_particles} min={min_particles} \
    Bias=WLTM min_sweeps={min_sweeps} min_flatness={min_flatness} collect_flatness={collect_flatness} min_collect_sweeps=1
TrialAddRemove weight=2 particle_type=pt1
TrialAddRemove weight=2 particle_type=pt2
Log [write].csv
Tune [write]_tune.csv multistate=true stop_after_cycle=1
#To print xyz for each macrostate in separate files, add the following arguments to the "Movie" lines below: multistate true multistate_aggregate false
Movie [write]_eq.xyz stop_after_cycle=1
Movie [write].xyz start_after_cycle=1
Energy [write]_en.csv multistate=true start_after_cycle=1
ProfileCPU [write]_profile.csv
NumParticles [write]_num0.csv particle_type=pt1 multistate=true start_after_cycle=1
#HeatCapacity [write]_cv.csv multistate=true start_after_cycle=1
CriteriaWriter [write]_crit.csv
ExtensiveMoments [write]_ext.csv
CriteriaUpdater trials_per_update=1e5
Run until=complete

# continue until all simulations on the node are complete
WriteFileAndCheck sim={sim} sim_start={sim_start} sim_end={sim_end} file_prefix={prefix}n{node}s file_suffix=_finished.txt output_file={prefix}n{node}_terminate.txt
Run until_file_exists={prefix}n{node}_terminate.txt trials_per_file_check={tpc}
""".format(**params))

def post_process(params):
    """ Skip the following checks if temperature is not the default """
    if np.abs(params['beta'] - 1./258.15) > 1e-5:
        return
    lnpi=macrostate_distribution.splice_files(prefix=params['prefix']+'n0s', suffix='_crit.csv', shift=False)
    #lnpi.reweight(delta_beta_mu=-4.5, inplace=True)
    num0 = multistate_accumulator.splice(prefix=params['prefix']+'n0s', suffix='_num0.csv')
    #num0.to_csv('num0.csv')
    #plt.plot(num0['average'])
    #plt.show()
    #print('num0', num0)
    lnpi.concat_dataframe(dataframe=num0, add_prefix='num0_')
    #lnpi.reweight(delta_beta_mu=-4.5, inplace=True)
    #lnpi.plot(show=True) # show before reweighting to equilibrium
    lnpi.equilibrium(delta_beta_mu_guess=-4.5)
    #lnpi.plot(show=True)
    for index, phase in enumerate(lnpi.split()):
        print('##############\nphase/index', index, '\n##############')
        n_gce = phase.average_macrostate()
        print('n_gce', n_gce)
        num0 = phase.ensemble_average('num0_average')
        print('num0', num0)
        if index == 0:
            assert np.abs(n_gce - 27) < 4
            pressure = -phase.ln_prob()[0]/params['beta']/params['cubic_side_length']**3
            print('pressure(K/Ang^3)', pressure)
            # convert pressure in K/Ang^3 to MPa
            # K/Ang^3 (1e30 A^3 / m^3)(Pa m^3/J)(MPa/1e6/Pa)(kB J/K)
            pres_conv = 1e30/1e6*physical_constants.BoltzmannConstant().value()
            print('pressure(MPa)', pressure*pres_conv)
            assert np.abs(pressure*pres_conv - 5.38) < 0.4
            assert np.abs(num0 - 14.8) < 4
        else:
            assert np.abs(n_gce - 183.8) < 4
            assert np.abs(num0 - 172.8) < 4
        #plt.plot(phase.dataframe()['num0_average'])
        #plt.show()
        print('mole frac0', num0/n_gce)
    plt.plot(lnpi.dataframe()['state'], lnpi.dataframe()['ln_prob'], label='FEASST')
    plt.xlabel('number of particles', fontsize=16)
    plt.ylabel('ln probability', fontsize=16)
    plt.legend(fontsize=16)
    #plt.savefig(params['prefix']+'_lnpi.png', bbox_inches='tight', transparent='True')
    #plt.show()

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
