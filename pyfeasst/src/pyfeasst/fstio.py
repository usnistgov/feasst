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
    params['queue_command'] = "sbatch --array={queue_task}-{max_restarts}%1 {prefix}_slurm.txt".format(**params)
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
                """echo "Server port {server_port} buffer_size {buffer_size}" |""".format(**params) + args.feasst_install+"""bin/fst > {prefix}{sim:03d}_fstout.txt""".format(**params),
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
    file_name = "{}{:03d}".format(params['prefix'],sim)
    write_feasst_script(params, script_file=file_name+'_fstin.txt')
    syscode = subprocess.call(
        args.feasst_install+'bin/fst < '+file_name+'_fstin.txt  > '+file_name+'_fstout.txt',
        shell=True, executable='/bin/bash')
    return syscode

def run_single(sim, params, args, sim_node_dependent_params, write_feasst_script, post_process):
    """ Run a single simulation. If all simulations are complete, run PostProcess. """
    if args.queue_task == 0:
        params['sim'] = sim
        syscode = seed_param_write_run(sim, params, args, sim_node_dependent_params, write_feasst_script)
    else: # if queue_task < 1, restart from checkpoint
        syscode = subprocess.call("echo \"Restart {}{:03d}_checkpoint.fst\" | {}bin/fst".format(params['prefix'], sim, args.feasst_install), shell=True, executable='/bin/bash')
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

def num_sites_in_fstprt(fstprt_file, feasst_install=""):
    """
    Return the integer number of sites in a fstprt file

    :param str fstprt_file:
        The name of the fstprt file which describes a particle.
    :param str feasst_install:
        An option to provide the path to the feasst build directory in order to convert a
        fstprt_file which begins with the characters '/feasst' to the feasst directory.

    >>> from pyfeasst import fstio
    >>> fstio.num_sites_in_fstprt('../../../particle/lj.txt')
    1
    >>> fstio.num_sites_in_fstprt('../../../particle/n-butane.txt')
    4
    >>> fstio.num_sites_in_fstprt('/feasst/particle/n-butane.txt', '../../../build/')
    4
    """
    if fstprt_file[:7] == '/feasst':
        if feasst_install == "":
            print('feasst_install:', feasst_install, 'should be used given as the path to the',
                  'feasst build directory in order to convert /feasst in the fstprt_file')
        fstprt_file = feasst_install + '..' + fstprt_file[7:]
    with open(fstprt_file, 'r') as file1:
        lines = file1.read().splitlines()
    num_sites = 0
    phase = 0
    for line in lines:
        #print('line', line)
        values = line.split()
        #print('values', values)
        if len(values) == 0 and phase == 2:
            phase = 3
        if phase == 2:
            num_sites += 1
        if len(values) == 0 and phase == 1:
            phase = 2
        if len(values) == 1:
            if values[0] == "Sites":
                phase = 1
    if num_sites == 0:
        print('reached the end of file without finding number of sites in', fstprt_file)
        assert False
    return num_sites

# For use in function below only
def write_bond_(btype, file1, descript, num_steps_override=-1):
    import copy
    if num_steps_override > 0:
        descript = copy.deepcopy(descript)
        if num_steps_override > 1:
            descript['num_steps'] = ' num_steps='+str(num_steps_override)
        else:
            descript['num_steps'] = ''
    if btype == 'bond':
        file1.write("""bond=true mobile_site={site0} anchor_site={site1}{num_steps}{reference_index}\n""".format(**descript))
    elif btype == 'angle':
        file1.write("""angle=true mobile_site={site0} anchor_site={site1} anchor_site2={site2}{num_steps}{reference_index}\n""".format(**descript))
    elif btype == 'dihedral':
        file1.write("""dihedral=true mobile_site={site0} anchor_site={site1} anchor_site2={site2} anchor_site3={site3}{num_steps}{reference_index}\n""".format(**descript))
    else:
        print('unrecognized bond type', btype)
    return True

def write_linear_grow_file(filename, num_sites=None, gce=0, reference_index=0, num_steps=4, base_weight=1, conf=0, conf2=-1, particle_type=0, angle=True, dihedral=True, particle_file=None, feasst_install=None):
    """
    Write TrialGrowFile input files for linear chain particles.

    :param str filename:
        The name of the file to write the TrialGrowFile text.
    :param int num_sites:
        The number of sites in the linear chain particle.
        If None, obtain num_sites from the particle_file.
    :param int gce:
        If 0, write only canonical ensemble trials.
        If 1, write only grand canonical ensemble insertion and deletion trials.
        If 2, write only grand canonical ensemble insertion trials
        If 3, wirte only Gibbs particle transfer trials.
    :param int reference_index:
        The index of the reference potential to be used for dual-cut configurational bias.
        If -1, do not use a reference potential.
    :param int num_steps:
        The number of steps to be used for (dual-cut) configurational bias.
    :param float base_weight:
        The weight used for determining the probability of attempting a trial.
    :param int conf:
        The configuration index for the trial (see System).
    :param int conf2:
        The second configuration index for the trial used with gibbs transfers only.
    :param int particle_type:
        Type of particle in configuration.
    :param bool angle:
        True if there are angle potentials present in the linear chain.
    :param bool dihedral:
        True if there are dihedral potentials present in the linear chain.
    :param str particle_file:
        If not None, use the particle_file to find the site names.
    :param str feasst_install:
        An option to provide the path to the feasst build directory in order to convert a
        particle_file which begins with the characters '/feasst' to the feasst directory.

    The weight for trial moves is the base weight, but reptations are divided by 4 and partial regrowth
    weights are divided by the number of stages.

    Reptations has been disabled because there is an issue with them.

    >>> from pyfeasst import fstio
    >>> fstio.write_linear_grow_file("grow_grand_canonical.txt", num_sites=3, gce=1, reference_index=0, num_steps=4)
    >>> with open('grow_grand_canonical.txt') as file1:
    ...     lines = file1.readlines()
    >>> lines[0]
    'TrialGrowFile\\n'
    >>> lines[1]
    '\\n'
    >>> lines[2]
    'particle_type=0 weight=1 transfer=true site=2 num_steps=4 reference_index=0\\n'
    >>> fstio.write_linear_grow_file("grow_grand_canonical.txt", gce=1, reference_index=0, num_steps=4, particle_file='../../../particle/spce_new.txt')
    >>> with open('grow_grand_canonical.txt') as file1:
    ...     lines = file1.readlines()
    >>> lines[2]
    'particle_type=0 weight=1 transfer=true site=H2 num_steps=4 reference_index=0\\n'
    """
    descript = {'reference_index':'', 'num_steps':'', 'weight':base_weight, 'conf':'',
                'particle_type':particle_type, 'conf2':conf2, 'rept_weight':base_weight/4.}
    if reference_index >= 0:
        descript['reference_index'] = ' reference_index='+str(reference_index)
    if num_steps > 1:
        descript['num_steps'] = ' num_steps='+str(num_steps)
    if conf > 0:
        descript['conf'] = ' configuration_index='+str(conf)
    site_names = list()
    if particle_file != None:
        if num_sites == None:
            num_sites = num_sites_in_fstprt(particle_file, feasst_install)
        # obtain site_names from particle_file
        if particle_file[:7] == '/feasst':
            if feasst_install == "":
                print('feasst_install:', feasst_install, 'should be used given as the path to the',
                      'feasst build directory in order to convert /feasst in the particle_file')
            particle_file = feasst_install + '..' + particle_file[7:]
        with open(particle_file, 'r') as file1:
            lines = file1.readlines()
        iline = 0
        while lines[iline] != 'Sites\n':
            iline += 1
            if iline > 1e6 or iline == len(lines):
                print('Cannot find "Sites" in particle_file', particle_file)
                assert False
        for isite in range(num_sites):
            site_names.append(lines[iline+2+isite].split(' ')[0])
    else:
        site_names = list(range(num_sites))
    #print('site_names', site_names)
    with open(filename, 'w') as file1:
        file1.write("TrialGrowFile\n\n")
        for inv in [True, False]:
            for trial_type in range(3+int(num_sites/2)): # 0: reptate, 1: full regrow, 2+: partial regrow
                wrote = False
                for site in range(num_sites):
                    for i in range(4): # build up to dihedrals
                        sign = -1
                        if trial_type == 0 and site != num_sites - 1:
                            sign = 1
                        if inv:
                            index = num_sites - site - 1 - sign*i
                        else:
                            index = site + sign*i
                        if index >= 0 and index < len(site_names):
                            #print('index', index, 'site_names', site_names)
                            descript['site'+str(i)] = site_names[index]

                    # full regrowth insertion/deletion
                    if trial_type == 1 and gce > 0:
                        if site == 0:
                            if gce == 3:
                                file1.write("""particle_type={particle_type}{conf} configuration_index2={conf2} weight={weight} gibbs_transfer=true site={site0}{num_steps}{reference_index} print_num_accepted=true\n""".format(**descript))
                            elif gce == 2:
                                file1.write("""particle_type={particle_type} weight={weight} add=true site={site0}{num_steps}{reference_index}{conf}\n""".format(**descript))
                            elif gce == 1:
                                file1.write("""particle_type={particle_type} weight={weight} transfer=true site={site0}{num_steps}{reference_index}{conf}\n""".format(**descript))
                            wrote = True
                        elif site == 1:
                            wrote = write_bond_('bond', file1, descript)
                        elif site == 2:
                            if angle:
                                wrote = write_bond_('angle', file1, descript)
                            else:
                                wrote = write_bond_('bond', file1, descript)
                        else:
                            if dihedral:
                                wrote = write_bond_('dihedral', file1, descript)
                            else:
                                wrote = write_bond_('bond', file1, descript)

#                    # reptation (must have num_steps 1)
#                    elif trial_type == 0 and gce == 0:
#                        if site == num_sites - 1:
#                            if num_sites == 2:
#                                wrote = write_bond_('bond', file1, descript, num_steps_override=1)
#                            elif num_sites == 3:
#                                if angle:
#                                    wrote = write_bond_('angle', file1, descript, num_steps_override=1)
#                                else:
#                                    wrote = write_bond_('bond', file1, descript, num_steps_override=1)
#                            elif num_sites > 3:
#                                if dihedral:
#                                    wrote = write_bond_('dihedral', file1, descript, num_steps_override=1)
#                                else:
#                                    wrote = write_bond_('bond', file1, descript, num_steps_override=1)
#                            else:
#                                print('unrecognized num_sites', num_sites)
#                                assert False
#                        else:
#                            if site == 0:
#                                file1.write("""particle_type {particle_type} weight {rept_weight}{conf} """.format(**descript))
#                            file1.write("""reptate true mobile_site {site0} anchor_site {site1}{reference_index}\n""".format(**descript))
#                            wrote = True

                    # partial regrow
                    if gce == 0 and trial_type > 1:
                        num_grow = trial_type - 1
                        if num_sites - site < num_grow:
                            if num_sites - site == num_grow - 1:
                                file1.write('particle_type='+str(descript['particle_type'])+' weight='+str(descript['weight']/(trial_type-2))+descript['conf']+' ')
                                wrote = True
                            if site == 1:
                                wrote = write_bond_('bond', file1, descript)
                            elif site == 2:
                                if angle:
                                    wrote = write_bond_('angle', file1, descript)
                                else:
                                    wrote = write_bond_('bond', file1, descript)
                            elif site != 0:
                                if dihedral:
                                    wrote = write_bond_('dihedral', file1, descript)
                                else:
                                    wrote = write_bond_('bond', file1, descript)

                if wrote:
                    file1.write("\n")

if __name__ == "__main__":
    import doctest
    doctest.testmod()
