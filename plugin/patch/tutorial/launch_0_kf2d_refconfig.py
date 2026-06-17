"""
Example reference configuration a 2D Kern-Frenkel patchy particle.
https://doi.org/10.1063/1.1569473
"""

import argparse
import numpy as np
from feasst import fstio

def parse():
    """ Parse arguments from command line or change their default values. """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--fstprt', type=str, default='/feasst/plugin/patch/particle/one_patch2d.txt',help='FEASST particle definition')
    parser.add_argument('--cutoff', type=float, default=1.5, help='the square well attractive interaction cutoff distance between centers')
    parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
    parser.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
    parser.add_argument('--num_jobs', type=int, default=1, help='Number of jobs in queue')
    parser.add_argument('--procs_per_job', type=int, default=1, help='number of processors')
    parser.add_argument('--run_type', '-r', type=int, default=0,
                        help='0: run, 1: submit to queue, 2: post-process')
    parser.add_argument('--seed', type=int, default=-1,
                        help='Random number generator seed. If -1, assign random seed to each sim.')
    parser.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
    parser.add_argument('--scratch', type=str, default=None,
                        help='Optionally write scheduled job to scratch/logname/jobid.')
    parser.add_argument('--queue_flags', type=str, default="", help='extra flags for queue (e.g., for slurm, "-p queue")')
    parser.add_argument('--job', type=int, default=0, help='job ID')
    parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
    parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

    # Convert arguments into a parameter dictionary, and add argument-dependent parameters.
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    params = vars(args)
    params['script'] = __file__
    params['prefix'] = 'kf2dref'
    params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
    params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
    params['procs_per_sim'] = 1
    params['num_sims'] = params['num_jobs']*params['procs_per_job']
    params['polar'] = np.pi/2 - 1e-8
    write_xyz(params, params['prefix']+'_in.xyz')
    params['polar'] = np.pi/2 + 1e-8
    write_xyz(params, params['prefix']+'_out.xyz')
    return params, args

def write_xyz(params, output_file):
    #now make this work for various polars and test
    r = 1.25
    params['x'] = r*np.cos(params['polar'])
    params['y'] = r*np.sin(params['polar'])
    params['xp'] = params['x'] - 1.
    with open(output_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""4
-1 10 10
0 0 0
1 1 0
0 {x} {y}
1 {xp} {y}""".format(**params))

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed={seed}
Configuration xyz_file={prefix}_in.xyz \
    particle_type=fluid:{fstprt} add_num_fluid_particles=2 \
    group=centers centers_site_type=A \
    patch_angleB=90 cutoff={cutoff}
Potential Model=HardSphere group=centers
Potential Model=SquareWell VisitModelInner=VisitModelInnerPatch group=centers
ThermoParams beta=1
Metropolis
Log output_file={prefix}{sim:03d}.csv clear_file=true
Run num_trials=1

MonteCarlo
RandomMT19937 seed={seed}
Configuration xyz_file={prefix}_out.xyz \
    particle_type=fluid:{fstprt} add_num_fluid_particles=2 \
    group=centers centers_site_type=A \
    patch_angleB=90 cutoff={cutoff}
Potential Model=SquareWell VisitModelInner=VisitModelInnerPatch group=centers
ThermoParams beta=1
Metropolis
Log output_file={prefix}{sim:03d}.csv
Run num_trials=1
""".format(**params))

def post_process(params):
    import numpy as np
    import pandas as pd
    df = pd.read_csv(params['prefix']+'000.csv')
    #print(df)
    #print(df['energy'])
    assert np.abs(float(df['energy'][0]) + 1) < 1e-8
    assert np.abs(float(df['energy'][2])) < 1e-8

if __name__ == '__main__':
    parameters, arguments = parse()
    fstio.run_simulations(params=parameters,
                          sim_job_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_job,
                          args=arguments)
