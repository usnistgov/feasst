"""
Isothermal-isobaric ensemble Monte Carlo simulation of SPC/E particles.
"""

import argparse
import json
from pyfeasst import fstio
from pyfeasst import physical_constants

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',
                        help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--fstprt', type=str, default='/feasst/particle/spce.fstprt',
                        help='FEASST particle definition')
    parser.add_argument('--temperature', type=float, default=500, help='temperature in Kelvin')
    parser.add_argument('--num_particles', type=int, default=500, help='number of particles')
    #parser.add_argument('--pressures', type=json.loads, default='{"pressure":[5.790E+01]}',
    parser.add_argument('--pressures', type=json.loads, default='{"pressure":[5.790E+01,5.521E+02,1.277E+03,2.270E+03,3.290E+03,3.432E+03,3.586E+03,3.728E+03]}',
                        help='dictionary with a list of pressures in units of atm.')
    parser.add_argument('--initial_cubic_side_length', type=int, default=25, help='cubic periodic boundary length')
    parser.add_argument('--tpc', type=int, default=int(1e4), help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=int(2e1), help='number of cycles for equilibraiton')
    parser.add_argument('--production', type=int, default=int(2e1), help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
    parser.add_argument('--run_type', '-r', type=int, default=0, help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--seed', type=int, default=-1,
                        help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
    parser.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
    parser.add_argument('--scratch', type=str, default=None,
                        help='Optionally write scheduled job to scratch/logname/jobid.')
    parser.add_argument('--queue_flags', type=str, default="", help='extra flags for queue (e.g., for slurm, "-p queue")')
    parser.add_argument('--node', type=int, default=0, help='node ID')
    parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
    parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

    # Convert arguments into a parameter dictionary, and add argument-dependent parameters.
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    params = vars(args)
    params['script'] = __file__
    params['prefix'] = 'npt'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_node'] = len(params['pressures']['pressure'])
    params['procs_per_sim'] = 1
    params['alpha'] = 5.6/params['initial_cubic_side_length']
    params['beta'] = 1./(params['temperature']*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
    params['num_sims'] = params['num_nodes']*params['procs_per_node']
    return params, args

def sim_node_dependent_params(params):
    """ Set parameters that depent upon the sim or node here. """
    params['pressure'] = params['pressures']['pressure'][params['sim']]
    # convert atm to energy/volume units (kJ/mol/A^3)
    params['pressure'] *= 101325/1e3/1e30*physical_constants.AvogadroConstant().value()

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {initial_cubic_side_length} particle_type0 {fstprt}
Potential VisitModel Ewald alpha {alpha} kmax_squared 38
Potential Model ModelTwoBodyFactory model0 LennardJones model1 ChargeScreened VisitModel VisitModelCutoffOuter erfc_table_size 2e4
Potential Model ChargeScreenedIntra VisitModel VisitModelBond
Potential Model ChargeSelf
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential 10
Metropolis
TrialTranslate tunable_param 2
TrialParticlePivot weight 0.5 particle_type 0 tunable_param 0.5
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}
CheckEnergy trials_per_update {tpc} decimal_places 4

# gcmc initialization
TrialAdd particle_type 0
Log trials_per_write {tpc} output_file {prefix}{sim}_init.csv
Tune
Run until_num_particles {num_particles}
Remove name0 TrialAdd name1 Log name2 Tune

# npt equilibration
ThermoParams beta {beta} pressure {pressure}
Metropolis trials_per_cycle {tpc} cycles_to_complete {equilibration}
TrialVolume weight 0.05 tunable_param 0.2 tunable_target_acceptance 0.5
Tune trials_per_tune 20
Log     trials_per_write {tpc} output_file {prefix}{sim}_eq.csv
Movie   trials_per_write {tpc} output_file {prefix}{sim}_eq.xyz
Density trials_per_write {tpc} output_file {prefix}{sim}_density_eq.csv
Run until complete
Remove name0 Tune name1 Log name2 Movie name3 Density

# npt production
Metropolis trials_per_cycle {tpc} cycles_to_complete {production}
Log        trials_per_write {tpc} output_file {prefix}{sim}.csv
Movie      trials_per_write {tpc} output_file {prefix}{sim}.xyz start_after_cycle 1
Energy     trials_per_write {tpc} output_file {prefix}{sim}_en.csv
Density    trials_per_write {tpc} output_file {prefix}{sim}_density.csv
Volume     trials_per_write {tpc} output_file {prefix}{sim}_volume.csv
ProfileCPU trials_per_write {tpc} output_file {prefix}{sim}_cpu.csv
Run until complete
""".format(**params))

def post_process(params):
    """
    Compare with https://www.nist.gov/mml/csd/informatics/lammps-md-equation-state-pressure-vs-density-spce-water
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    ens = np.zeros(shape=(params['num_sims'], 2))
    rhos = np.zeros(shape=(params['num_sims'], 2))
    dens_conv = 18.01528*1e30/1e3/physical_constants.AvogadroConstant().value()
    for sim in range(params['num_sims']):
        log = pd.read_csv(params['prefix']+str(sim)+'.csv')
        #assert int(log['num_particles_of_type0'][0]) == params['num_particles']
        energy = pd.read_csv(params['prefix']+str(sim)+'_en.csv')
        ens[sim] = np.array([energy['average'][0],
                             energy['block_stdev'][0]])/params['num_particles']
        density = pd.read_csv(params['prefix']+str(sim)+'_density.csv')
        #print('density', density)
        rhos[sim] = np.array([density['average'][0],
                              density['block_stdev'][0]])
        rhos[sim] *= dens_conv # molec/A^3 to kg/m^3
        #print('rhos[sim]', rhos[sim])
        ens[sim] /= 4.184 # kJ/mol to kcal/mol
    print('rhos(kg/m^3)', rhos)
    print('ens(kcal/mol)', ens)
    assert np.abs(rhos[0][0] - 800) < 55
    assert np.abs(rhos[1][0] - 850) < 35
    assert np.abs(rhos[2][0] - 900) < 35
    assert np.abs(rhos[3][0] - 950) < 35
    assert np.abs(rhos[4][0] - 990) < 35
    assert np.abs(rhos[5][0] - 995) < 35
    assert np.abs(rhos[6][0] - 1000) < 35
    assert np.abs(rhos[7][0] - 1005) < 35
    assert np.abs(ens[0][0] + 8.362) < 0.45
    assert np.abs(ens[1][0] + 8.565) < 0.45
    assert np.abs(ens[2][0] + 8.758) < 0.45
    assert np.abs(ens[3][0] + 8.934) < 0.45
    assert np.abs(ens[4][0] + 9.064) < 0.45
    assert np.abs(ens[5][0] + 9.078) < 0.45
    assert np.abs(ens[6][0] + 9.095) < 0.45
    assert np.abs(ens[7][0] + 9.108) < 0.45

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_node_dependent_params=sim_node_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=arguments)
