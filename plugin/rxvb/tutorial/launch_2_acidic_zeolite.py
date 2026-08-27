"""
particle types:
(0) OH <---> (2) O- (fixed)
(1) H2O <---> (3) H3O (fluid)
(4) framework (fixed, unreactive)
(5) blocks

site types:
0: O in OH particle type 0 (AV pair 1)
1: H in OH particle type 0
2: O in H2O particle type 1 (AV pair 1)
3: H in H2O particle type 1
4: Dummy in H2O particle type 1
5: O in O- particle type 2 (AV pair 2)
6: Dummy in O- particle type 2
7: O in H3O particle type 3 (AV pair 2)
8: H in H3O particle type 3
9: framework
10: framework
11: framework

Start with 1 OH, num_particles H2O and framework
"""

import argparse
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from feasst import fstio
from feasst import macrostate_distribution
from feasst import physical_constants
from pathlib import Path


def parse():
    # Parse arguments from command line or change their default values.
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--feasst_exe', type=str, default='feasst', help='FEASST executable')
    parser.add_argument('--temperature', type=float, default=373, help='temperature (K)')
    parser.add_argument('--mu', type=float, default=-1000, help='chemical potential of fixed site HS (adsorbed)')
    parser.add_argument('--mu_init', type=float, default=10, help='initial chemical potential')
    parser.add_argument('--num_particles', type=int, default=1, help='minimum number of particles')
    #parser.add_argument('--num_particles', type=int, default=250, help='minimum number of particles')
    parser.add_argument('--weight', type=float, default=0.01, help='weight for rxn trials')
    parser.add_argument('--tpc', type=int, default=int(1e6), help='trials per cycle')
    parser.add_argument('--equilibration', type=int, default=1e4, help='number of cycles for equilibration')
    parser.add_argument('--production', type=int, default=1e6, help='number of cycles for production')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
    parser.add_argument('--procs_per_job', type=int, default=1, help='number of processors')
    parser.add_argument('--run_type', '-r', type=int, default=0,
                        help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--seed', type=int, default=-1,
                        help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--max_restarts', type=int, default=0, help='Number of restarts in queue')
    parser.add_argument('--num_jobs', type=int, default=1, help='Number of jobs in queue')
    parser.add_argument('--scratch', type=str, default=None,
                        help='Optionally write scheduled job to scratch/logname/jobid.')
    parser.add_argument('--job', type=int, default=0, help='job ID')
    parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
    parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')
    parser.add_argument('--out', type=str, default=None, help='Output directory for this run')

    # Convert arguments into a parameter dictionary, and add argument-dependent parameters.
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    params = vars(args)
    params['script'] = __file__
    params['prefix'] = 'acid_zeo'
    params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.95*params['hours_terminate'] - 0.05 # terminate FEASST before SLURM
    params['procs_per_sim'] = 1
    params['num_sims'] = params['procs_per_job']*params['num_jobs']
    params['beta'] = 1./(params['temperature']*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
    params['avb_outer'] = 3.7
    params['avb_inner'] = 2.4
    temp_tag = f"temp{int(args.temperature)}"
    nw_tag = f"{args.num_particles}w"
    mu_tag = f"mu_{abs(int(args.mu))}"

    if args.out:
        outdir = Path(args.out)
    else:
        outdir = Path('runs') / temp_tag / nw_tag / mu_tag

    outdir.mkdir(parents=True, exist_ok=True)
    params['outdir'] = str(outdir)

    params['assets_dir'] = str(Path(__file__).parent.resolve())

    base_prefix = params['prefix']                     # "acid_zeo"
    params['prefix'] = str(outdir / base_prefix)       # "runs/temp_/w/mu_xxx/acid_zeo"
    params['sim_id_file'] = params['prefix'] + '_sim_ids.txt'

    params['assets_dir'] = str(Path(__file__).parent.resolve())
    return params, args

def sim_job_dependent_params(params):
    """ Define parameters that are dependent on the sim or job. """
    params['sim_start'] = 0
    params['sim_end'] = params['num_sims'] - 1
    params['morph'] = params['sim'] % 2
    if params['morph'] == 1:
        params['rxn']="""TrialMorph weight={weight} particle_type=OH,H2O particle_type_morph=O,H3O print_num_accepted=true reference_index=0""".format(**params)
        params['two_body_potential'] = """Potential Model=ModelTwoBodyFactory models=LennardJonesForceShift,ChargeScreened erfc_table_size=2e4""".format(**params)
    else:
        params['rxn'] = """TrialRxVB weight={weight} target_particle_type=OH target_site=0 target_particle_type_morph=O particle_type=H2O site=2 particle_type_morph=H3O print_num_accepted=true reference_index=0""".format(**params)
        params['two_body_potential'] = """NeighborCriteria site_type0=0 site_type1=2 site_type0_alt=5 site_type1_alt=7 maximum_distance={avb_outer} minimum_distance={avb_inner} potential_index=1
Potential Model=ModelTwoBodyFactory models=LennardJonesForceShift,ChargeScreened erfc_table_size=2e4 EnergyMap=EnergyMapNeighborCriteria""".format(**params)

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed={seed}
Configuration side_length=20.022,19.899,13.383 \
    particle_type=OH:{assets_dir}/oh.fstprt,H2O:{assets_dir}/spce_with_dummy.fstprt,O:{assets_dir}/o-.fstprt,H3O:{assets_dir}/h3o.fstprt,framework:{assets_dir}/framework.fstprt,blocks:{assets_dir}/block.fstprt \
    group=mobile,framework,blocks mobile_particle_type=OH,H2O,O,H3O \
    framework_particle_type=framework \
    blocks_particle_type=blocks \
    add_num_OH_particles=1 add_num_framework_particles=1
Potential VisitModel=Ewald tolerance=1e-4 tolerance_num_sites=2000
{two_body_potential}
Potential Model=ChargeScreenedIntra VisitModel=VisitModelBond
Potential Model=ChargeSelf
Potential Model=ModelHardShape shape_file={assets_dir}/MFI_Zeolite_111_1-54_blocks.txt cavity=false group=blocks
ZeroBackground
RefPotential VisitModel=DontVisitModel
ThermoParams beta={beta} chemical_potential={mu_init},{mu_init}
Metropolis
TrialTranslate weight_per_number_fraction=1 number_fraction_exclude_type0=4 particle_type=H2O
TrialTranslate weight_per_number_fraction=1 number_fraction_exclude_type0=4 particle_type=H3O
TrialParticlePivot weight_per_number_fraction=1 number_fraction_exclude_type0=4 particle_type=OH
TrialParticlePivot weight_per_number_fraction=1 number_fraction_exclude_type0=4 particle_type=H2O
TrialParticlePivot weight_per_number_fraction=1 number_fraction_exclude_type0=4 particle_type=O
TrialParticlePivot weight_per_number_fraction=1 number_fraction_exclude_type0=4 particle_type=H3O
TrialSwapPosition weight=0.1 particle_types=H2O,H3O
CheckEnergy trials_per_update={tpc} decimal_places=6
Checkpoint checkpoint_file={prefix}{sim}_checkpoint.fst num_hours={hours_checkpoint} num_hours_terminate={hours_terminate}

# write framework once for visualization
Movie output_file={prefix}n{job}s{sim}_framework.xyz group=framework clear_file=true
WriteStepper analyze_name=Movie
Remove name=Movie

# gcmc initialization
Let [write]=trials_per_write={tpc} output_file {prefix}n{job}s{sim}
Log [write]_fill.csv
Tune
Run until_num_particles={num_particles} particle_type=H2O Trial=TrialAdd
WriteStepper analyze_name=Log
Remove name=Log

# equilibration
ThermoParams beta={beta} chemical_potential=0,0,0,{mu}
{rxn}
Metropolis trials_per_cycle={tpc} cycles_to_complete={equilibration}
Log [write]_eq.csv
Movie [write]_eq.xyz group=mobile
NumParticles [write]_eq_num2.csv append=true rewrite_header=false particle_type=O
NumParticles [write]_eq_num3.csv append=true rewrite_header=false particle_type=H3O
Run until=complete
Remove name=Tune,Log,Movie,NumParticles,NumParticles

# production
Metropolis trials_per_cycle={tpc} cycles_to_complete={production}
Log [write].csv
Movie [write].xyz group=mobile
Energy [write]_en.csv append=true rewrite_header=false
NumParticles [write]_num2.csv append=true rewrite_header=false particle_type=O
NumParticles [write]_num3.csv append=true rewrite_header=false particle_type=H3O
CPUTime [write]_cpu.csv append=true
ProfileCPU [write]_profile.csv
Run until=complete

# continue until all simulations on the job are complete
WriteFileAndCheck sim={sim} sim_start={sim_start} sim_end={sim_end} file_prefix={prefix}n{job}s file_suffix=_finished.txt output_file={prefix}n{job}_terminate.txt
Run until_file_exists={prefix}n{job}_terminate.txt trials_per_file_check={tpc}
""".format(**params))

def post_process(params):
    assert True # placeholder

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_job_dependent_params=sim_job_dependent_params,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_job,
                          args=arguments)
