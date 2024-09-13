"""
This module provides some utility input / output functions for use with the FEASST simulation program.
"""

import os
import time
import sys
import random
import json
#import socket
import subprocess
from multiprocessing import Pool
import multiprocessing
from itertools import repeat
import numpy as np
from pathlib import Path

def dict_to_argparse(dictionary):
    """
    Converts a dictionary to a string that argparse may read from command line

    >>> from pyfeasst import fstio
    >>> params = {'run_type': 0, 'feasst_install': '/path/to/feasst/', 'queue_flags': ""}
    >>> fstio.dict_to_argparse(params)
    ' --run_type 0 --feasst_install /path/to/feasst/ --queue_flags ""'
    """
    rtrn = str()
    for label in dictionary:
        if dictionary[label] == "":
            rtrn += ' --'+label+' ""'
        else:
            rtrn += ' --'+label+' '+str(dictionary[label])
    return rtrn

def vector3d_to_list(vec):
    """
    Converts a swig stl vector to python list

    >>> from pyfeasst import fstio
    >>> fstio.vector3d_to_list([[[[0]]]])
    [[[[0]]]]
    """
    lst = list()
    for _, vec1 in enumerate(vec):
        lst2 = list()
        for _, vec2 in enumerate(vec1):
            lst3 = list()
            for _, vec3 in enumerate(vec2):
                lst3.append(vec3)
            lst2.append(lst3)
        lst.append(lst2)
    return lst

def read_checkpoint(filename):
    """
    Return contents of checkpoint file as a string

    >>> from pyfeasst import fstio
    >>> table = fstio.read_checkpoint('../../tests/tutorial_0_table.txt')
    >>> table[:13]
    '6867 11 11 11'
    """
    with open (filename, "r") as myfile:
        checkpoint=myfile.readlines()
    assert(len(checkpoint) == 1)  # checkpoint files should have only one line
    return checkpoint[0]

def all_sims_complete(filename, num_sims):
    """
    Read filename and see if all sim ID's from [0, num_sims-1] are present (e.g., complete).
    If no file is present, also consider the simulation incomplete.

    >>> from pyfeasst import fstio
    >>> all_sims_complete('../../tests/lj_sim_ids.txt', 8)
    True
    >>> all_sims_complete('../../tests/lj_sim_ids2.txt', 8)
    False
    >>> all_sims_complete('../../tests/not_a_file.txt', 8)
    False
    """
    if not os.path.isfile(filename):
        return False
    try:
        with open(filename, 'r') as file1:
            lines = file1.read().splitlines()
    except FileNotFoundError:
        return False
    ids = map(str, list(range(num_sims)))
    for line in lines:
        ids = [i for i in ids if i not in line]
    if len(ids) == 0:
        return True
    return False

def slurm_single_node(params):
    """
    Write slurm script to fill one node.

    :param dict params:
        Must have the following keys:
        procs_per_node: number of processors per node,
        procs_per_sim: number of processors per simulation,
        minutes: maximum number of minutes for job in queue,
        prefix: prefix for all output file names,
        script: script file name,
        sim_id_file: filename to write simulation id's for later checking of status,
        max_restarts: maximum number of restarts,
        node: node index.
        Optional keys:
        scratch: location of local scratch space on HPC node (disabled by default),
        scratch_hours_per_sync: hours between synching scratch to destination (default: 5).

    This function also adds the key 'queue_command' to the params dictionary,
    which is assumed to output the job id.
    """
    if not 'max_restarts' in params:
        params['max_restarts'] = '0'
    params['queue_command'] = "sbatch --array=0-" + str(params['max_restarts']) + "%1 " + params['prefix'] + "_slurm.txt"
    params['command_to_queue_id'] = " | tail -1 | awk '{print $4}'"
    if params['scratch'] == None or params['scratch'] == 'None':
        params['scratch_slurm_preamble'] = ''
        params['scratch_slurm_postamble'] = ''
    else:
        if 'scratch_hours_per_sync' not in params:
            params['scratch_seconds_per_sync'] = 5*60*60
        else:
            params['scratch_seconds_per_sync'] = 60*60*params['scratch_hours_per_sync']
        params['scratch_slurm_preamble'] = """original_dir=$PWD; echo $original_dir
scratch={scratch}/$LOGNAME/${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}/; mkdir -p $scratch; cd $scratch; echo "scratch:$scratch"
rsync -au $original_dir/* .
echo "while [ 1 -le 2 ]; do sleep {scratch_seconds_per_sync}; rsync -au * $original_dir/; done" > {prefix}_sync.sh
chmod a+x {prefix}_sync.sh
./{prefix}_sync.sh &
ls""".format(**params)
        params['scratch_slurm_postamble'] = 'rsync -au . $original_dir/'
    if not 'queue_flags' in params:
        params['queue_flags'] = ''
    with open(params['prefix'] + '_slurm.txt', 'w', encoding='utf-8') as myfile:
        myfile.write("""#!/bin/bash
#SBATCH -n {procs_per_node}
#SBATCH -N 1
#SBATCH -t {minutes}:00
#SBATCH -o {prefix}_slurm_%A_%a.txt
#SBATCH -e {prefix}_slurm_%A_%a.txt
#SBATCH {queue_flags}
echo "Running ID ${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}} on $(hostname) at $(date) in $PWD"
cd $PWD
{scratch_slurm_preamble}
export OMP_NUM_THREADS={procs_per_sim}
python {script} --run_type 0 --node {node} --queue_id $SLURM_ARRAY_JOB_ID --queue_task $SLURM_ARRAY_TASK_ID
{scratch_slurm_postamble}
if [ $? == 0 ] && [ ! -f {sim_id_file} ]; then
  echo "Simulation complete. Cancelling restarts."
  scancel $SLURM_ARRAY_JOB_ID
else
  echo "Job is terminating, to be restarted again"
fi
echo "Time is $(date)"
""".format(**params))

def server(params, args, sim_node_dependent_params, write_feasst_script):
    if 'port' not in params:
        params['port'] = 54321
    if 'buffer_size' not in params:
        params['buffer_size'] = 1000
    #print('starting server localhost:', params['port'])
    params['server_port'] = params['port'] + params['sim']
    if write_feasst_script == None:
        syscode = subprocess.call(
            """echo "Server port {server_port} buffer_size {buffer_size}" |""".format(**params) + args.feasst_install+'bin/fst > '+params['prefix']+str(params['sim'])+'_run.log',
            shell=True, executable='/bin/bash')
    else:
        syscode = seed_param_write_run(params['sim'], params, args, sim_node_dependent_params, write_feasst_script)
    return syscode

def reseed(params):
    #print('reseeding?', 'seed' in params)
    if 'seed' in params:
        if params['seed'] == -1:
            params['seed'] = random.randrange(int(1e9))
            #print('reseeding', params['seed'])

def seed_param_write_run(sim, params, args, sim_node_dependent_params, write_feasst_script):
    reseed(params)
    if sim_node_dependent_params:
        sim_node_dependent_params(params)
    file_name = params['prefix']+str(sim)+'_run'
    write_feasst_script(params, script_file=file_name+'.txt')
    syscode = subprocess.call(
        args.feasst_install+'bin/fst < '+file_name+'.txt  > '+file_name+'.log',
        shell=True, executable='/bin/bash')
    return syscode

def run_single(sim, params, args, sim_node_dependent_params, write_feasst_script, post_process):
    """ Run a single simulation. If all simulations are complete, run PostProcess. """
    if args.queue_task == 0:
        params['sim'] = sim
        syscode = seed_param_write_run(sim, params, args, sim_node_dependent_params, write_feasst_script)
    else: # if queue_task < 1, restart from checkpoint
        syscode = subprocess.call(
            args.feasst_install+'bin/rst --checkpoint_file '+params['prefix']+str(sim)+'_checkpoint.fst',
            shell=True, executable='/bin/bash')
    if syscode == 0: # if simulation finishes with no errors, write to sim id file
        with open(params['sim_id_file'], 'a', encoding='utf-8') as file1:
            file1.write(str(sim)+'\n')
        # if all sims are complete, post process or test once (by removing sim id file)
        if all_sims_complete(params['sim_id_file'], params['num_sims']):
            removed = False
            try:
                os.remove(params['sim_id_file'])
                removed = True
            except FileNotFoundError:
                pass
            if post_process is not None and removed:
                post_process(params)
    return syscode

def split(a, n):
    """
    From https://stackoverflow.com/questions/2130016/splitting-a-list-into-n-parts-of-approximately-equal-length

    >>> from pyfeasst import fstio
    >>> list(fstio.split(range(11), 3))
    [range(0, 4), range(4, 8), range(8, 11)]
    """
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

def read_write_params(args, params):
    #print('queue id', args.queue_id)
    if args.queue_id != -1: # if run from queue
        #print('args.queue_task', args.queue_task)
        if args.queue_task == 0: # read param file if not checkpoint
            with open(params['prefix']+'_params'+str(args.queue_id)+'.json', 'r') as file1:
                params = json.load(file1)
            #print('params', params)
    else:
        with open(params['prefix']+'_params.json', 'w') as file1:
            file1.write(json.dumps(params, indent=2))
    return params

def run_simulations(params, queue_function, args, write_feasst_script=None, client=None, sim_node_dependent_params=None, post_process=None):
    """
    Run a simulation either locally in the shell or queue on HPC nodes

    :param dict params:
        Must have the following keys:
        sim_id_file: filename to write simulation id's for later checking of status,
        prefix: prefix for all output file names,
        max_restarts: maximum number of restarts,
        procs_per_node: number of processors per node,
        procs_per_sim: number of processors per sim,
        num_nodes: number of nodes,
        node: node index.
    :param function sim_node_dependent_params:
        The name of the function that assigns parameters based on the sim and node.
        The only argument is the parameters.
    :param function write_feasst_script:
        The name of the function to write the feasst text interface file,
        which has the first argment as the parameters and the second argument as the filename.
    :param function client:
        If write_feasst_script is None, instead run in server mode with default port 54321 and default buffer_size 1000.
    :param function post_process:
        The name of the function to post process all simulations once complete,
        and has the only argument as the params.
    :param function queue_function:
        The name of the function to queue one node and has the only argument as the params.
    :param namespace args:
        Arguments from argparse.
    """
    with open(params['sim_id_file'], 'w') as file1:
        file1.close() # clear file, then append sim id when complete
    if args.run_type == 0: # run directly
        #print('getting params')
        params = read_write_params(args, params)
        #print('have params', params['beta'])
        if client is None:
            #print('Run text interface')
            with Pool(params['procs_per_node']) as pool:
                sims = list(split(range(params['num_sims']), params['num_nodes']))[params['node']]
                codes = pool.starmap(run_single, zip(sims, repeat(params), repeat(args),
                                                     repeat(sim_node_dependent_params),
                                                     repeat(write_feasst_script),
                                                     repeat(post_process)))
                if np.count_nonzero(codes) > 0:
                    sys.exit(1)
        else:
            #print('Run in server client mode')
            sims = list(split(range(params['num_sims']), params['num_nodes']))[params['node']]
            #procs = list()
            seed = 0
            #print('beta', params['beta'])
            if 'seed' in params:
                seed = params['seed']
            server_procs = list()
            client_procs = list()
            for s in sims:
                params['sim'] = s
                #print('s', s)
                proc = multiprocessing.Process(target=server, args=(params, args, sim_node_dependent_params, write_feasst_script))
                proc.start()
                server_procs.append(proc)
            time.sleep(1) # wait for servers to start before connecting client
            for s in sims:
                params['sim'] = s
                #print('s', s)
                if seed == -1:
                    params['seed'] = random.randrange(int(1e9))
                proc = multiprocessing.Process(target=client, args=(params,))
                proc.start()
                client_procs.append(proc)
            for proc in client_procs:
                proc.join()
            for proc in server_procs:
                proc.terminate()
            if post_process is not None:
                post_process(params)

    elif args.run_type == 1: # queue
        queue_id_file = params['prefix']+ '_queue_ids.txt'
        with open(queue_id_file, 'w') as file1:
            file1.close() # empty file contents
        for node in range(params['num_nodes']):
            params['node'] = node
            queue_function(params)
            subprocess.check_call(params['queue_command'] + params['command_to_queue_id'] + " >> " + queue_id_file, shell=True, executable='/bin/bash')
            with open(queue_id_file, 'r') as file1:
                queue_id = file1.read().splitlines()[-1]
            with open(params['prefix']+'_params'+queue_id+'.json', 'w') as file1:
                file1.write(json.dumps(params, indent=2))
    elif args.run_type == 2: # post process
        post_process(params)
    else:
        assert False  # unrecognized run_type

def combine_tables_two_rigid_body(prefix, suffix, num_procs, num_header_lines=6):
    """ Combine table files generated in parallel by TableTwoRigidBody3D."""
    lines = list()
    for proc in range(num_procs):
        filename = prefix + str(proc) + suffix
        with open(filename, 'r') as file1:
            lines.append(file1.readlines())
    with open(prefix + suffix, 'w') as file1:
        # write header
        for line in range(6):
            file1.write(lines[0][line])
        proc = 0
        done = 0
        line = 6
        itr = 0
        while done == 0:
            if line >= len(lines[proc]):
                done = 1
            else:
                file1.write(lines[proc][line])
            #print('proc', proc, 'line', line)
            proc += 1
            if proc == num_procs:
                proc = 0
                line += 1
            assert itr < 1e50

#def connect_to_port(params, max_attempts=5):
#    """Attempt to connect to port for server-client interface."""
#    port = params['port']+params['sim']
#    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#    connected = False
#    attempts = 0
#    max_attempts = 5
#    while not connected:
#        try:
#            sock.connect(("localhost", port))
#            connected = True
#        except:
#            time.sleep(1)
#            print('sim', params['sim'], 'failed to connect to port', port)
#            attempts += 1
#            assert attempts < max_attempts, "Could not connect to port."
#    return sock

if __name__ == "__main__":
    import doctest
    doctest.testmod()
