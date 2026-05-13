"""
Example ternary Gibbs ensemble Monte Carlo simulation of CO2 - methane - ethane with SAFT-based CG models.
"""

import argparse
from pyfeasst import fstio
from pyfeasst import physical_constants

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_install', type=str, default='../../../build/',help='FEASST install directory (e.g., the path to build)')
    parser.add_argument('--pt1', type=str, default='/feasst/particle/mie_C1H4_methane.txt', help='FEASST particle definition')
    parser.add_argument('--pt2', type=str, default='/feasst/particle/dimer_mie_CO2.txt', help='FEASST particle definition')
    parser.add_argument('--pt3', type=str, default='/feasst/particle/mie_C2H6_ethane.txt', help='FEASST particle definition')
    parser.add_argument('--beta', type=float, default=1./250.0, help='inverse temperature (K)')
    parser.add_argument('--pressure', type=float, default=3.03975, help='pressure (MPa)')
    parser.add_argument('--tpc', type=int, default=int(1e4), help='trials per cycle')
    parser.add_argument('--equil_npt', type=int, default=int(1e1), help='number of cycles for np equilibraiton')
    parser.add_argument('--equil', type=int, default=int(1e2), help='number of cycles for Gibbs equilibraiton')
    parser.add_argument('--production_cycles', type=int, default=int(1e2),help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=56, help='hours until termination')
    parser.add_argument('--procs_per_job', type=int, default=1, help='number of processors')
    parser.add_argument('--run_type', '-r', type=int, default=0,help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--seed', type=int, default=-1,help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
    parser.add_argument('--num_jobs', type=int, default=1, help='Number of jobs in queue')
    parser.add_argument('--scratch', type=str, default=None,help='Optionally write scheduled job to scratch/logname/jobid.')
    parser.add_argument('--queue_flags', type=str, default="", help='extra flags for queue (e.g., for slurm, "-p queue")')
    parser.add_argument('--job', type=int, default=0, help='job ID')
    parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
    parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

    # Convert arguments into a parameter dictionary, and add argument-dependent parameters.
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    params = vars(args)
    params['script'] = __file__
    params['prefix'] = 'ternary'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_jobs']*params['procs_per_job']
    # convert pressure in MPa to K/Ang^3
    # MPa(1e6 Pa/MPa)(J/m^3/Pa)(m^3/1e30/A^3)(K/kB/J)
    params['pressure'] = params['pressure']/1e24/physical_constants.BoltzmannConstant().value()
    params['dccb_cut'] = 3.5 # cutoff for dual-cut cb in liquid
    params['num_dccb'] = 4   # number of cb steps per site in liquid
    return params, args

def sim_job_dependent_params(params):
    """ Set parameters that depent upon the sim or job here. """
    for ptype, fstprt in enumerate([params['pt1'], params['pt2'], params['pt3']]):
        ptype_name = 'pt'+str(ptype+1)
        prepend = "{}{:03d}_p{}".format(params['prefix'], params['sim'], ptype+1)
        fstio.write_linear_grow_file(prepend+"_vapor_grow_add.txt", particle_type=ptype_name, gce=2, conf="vapor", ref="noixn", num_steps=1, particle_file=fstprt, feasst_install=params['feasst_install'])
        fstio.write_linear_grow_file(prepend+"_liquid_grow_add.txt", particle_type=ptype_name, gce=2, conf="liquid", ref="dccb", num_steps=params['num_dccb'], particle_file=fstprt, feasst_install=params['feasst_install'])
        fstio.write_linear_grow_file(prepend+"_grow_gibbs.txt", particle_type=ptype_name, gce=3, conf="vapor", conf2="liquid", ref="dccb", num_steps=params['num_dccb'], particle_file=fstprt, feasst_install=params['feasst_install'])

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed={seed}
For [config]:[length]=vapor:81,liquid:55
    Configuration name=[config] cubic_side_length=[length] particle_type=pt1:{pt1},pt2:{pt2},pt3:{pt3} \
        model_param_file=/feasst/particle/mie_model_parameters.txt
    Potential Model=Mie table_size=1e4 config=[config]
    Potential VisitModel=LongRangeCorrections config=[config]
    RefPotential ref=noixn VisitModel=DontVisitModel config=[config]
    WriteModelParams output_file={prefix}{sim:03d}_[config]_model_params.txt config=[config]
EndFor
RefPotential ref=dccb config=vapor VisitModel=DontVisitModel
RefPotential ref=dccb config=liquid Model=Mie table_size=1e4 VisitModel=VisitModelCell cutoff={dccb_cut} min_length={dccb_cut}
ThermoParams beta {beta} chemical_potential=5,5,5 pressure {pressure}
Metropolis
For [config]=vapor,liquid
    For [pt]=1,2,3
        TrialTranslate weight_per_number_fraction=0.2 particle_type=pt[pt] tunable_param=0.2 config=[config]
    EndFor
    For [pt]=2,3
        TrialParticlePivot particle_type=pt[pt] tunable_param=0.5 config=[config]
    EndFor
EndFor
Checkpoint checkpoint_file={prefix}{sim:03d}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# fill both boxes with particles
Let [write]=trials_per_write={tpc} output_file={prefix}{sim:03d}
Log [write]_fill.csv
Tune
For [config]:[num1]:[num2]:[num3]=vapor:330:330:410,liquid:100:450:700
    Movie [write]_[config]_fill.xyz config=[config]
    For [pt]:[num]=1:[num1],2:[num2],3:[num3]
        TrialGrowFile grow_file={prefix}{sim:03d}_p[pt]_[config]_grow_add.txt
        Run until_num_particles=[num] particle_type=pt[pt] config=[config]
        Remove name_contains=add
    EndFor
    Remove name=Movie
EndFor
Remove name=Tune,Log

# npt equilibrate both boxes
Metropolis trials_per_cycle={tpc} cycles_to_complete={equil_npt}
For [config]=vapor,liquid
    TrialVolume weight=0.005 tunable_param=0.2 tunable_target_acceptance=0.5 ref=noixn config=[config]
    Movie [write]_[config]_npt.xyz config=[config]
EndFor
Log [write]_npt.csv
ProfileCPU [write]_npt_profile.csv
CheckEnergy trials_per_update={tpc} decimal_places=8
# a new tune is required when new Trials are introduced
# decrease trials per due to infrequency of volume transfer attempts
Tune trials_per_tune=20
Run until=complete
Remove name=Tune,Log,Movie,Movie,ProfileCPU

# Gibbs equilibration
Metropolis trials_per_cycle={tpc} cycles_to_complete={equil}
For [pt]=1,2,3
    TrialGrowFile grow_file={prefix}{sim:03d}_p[pt]_grow_gibbs.txt
EndFor
Log [write]_eq.csv
For [config]=vapor,liquid
    Movie [write]_[config]_eq.xyz config=[config]
EndFor
ProfileCPU [write]_eq_profile.csv
Run until=complete
Remove name=Log,Movie,Movie,ProfileCPU

# Gibbs ensemble production
Metropolis trials_per_cycle={tpc} cycles_to_complete={production_cycles}
Log [write].csv
For [config]=vapor,liquid
    For [analyze]:[file]=Density:_dens.csv,Movie:.xyz,Energy:_en.csv,Volume:_vol.csv
        [analyze] [write]_[config][file] config=[config]
    EndFor
    For [pt]=1,2,3
        NumParticles [write]_[config]_n[pt].csv particle_type=pt[pt] config=[config]
    EndFor
EndFor
For [analyze]:[file]=ProfileCPU:_profile.csv,CPUTime:_cpu.csv
    [analyze] [write][file]
EndFor
Run until=complete
""".format(**params))

def post_process(params):
    import numpy as np
    import pandas as pd
    ns = list()
    nsa = list()
    sim = 0
    nexp = [267, 237, 244]
    frac = nexp/np.sum(nexp)
    for i,n in enumerate([1,2,3]):
        ns.append(pd.read_csv("{}{:03d}_vapor_n{}.csv".format(params['prefix'], sim, n)))
        nsa.append(ns[-1]['average'][0])
    for i,n in enumerate([1,2,3]):
        diff = np.abs(ns[i]['average'][0]/np.sum(nsa) - frac[i])
        #print('n', ns[-1])
        #print('diff', diff)
        if diff > 0.03:
            print('n', n, 'nexp', nexp[i], 'n', ns[-1], 'diff', diff)
            assert False

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_job_dependent_params=sim_job_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_job,
                          args=arguments)
