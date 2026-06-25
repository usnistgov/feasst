"""
Obtain information from FEASST scripts, particle files, etc
"""

import copy
from importlib.resources import files
#import logging as log

# Configure log to a file
#log.basicConfig(filename="file_reader.debug.txt", level=log.DEBUG, format='[%(filename)s:%(lineno)d] %(message)s')

def find_between(text: str, start: str, end: str):
    """
    Find the string in text that lies between the start and end strings.
    If the start is not found, return an empty string.
    If the start is found, but the end is not found, return everything from the start.

    >>> from feasst import file_reader
    >>> file_reader.find_between("The [deer] jumps over the fence.", '[', ']')
    'deer'
    >>> file_reader.find_between("The .deer. jumps over the fence.", ' .', ". ")
    'deer'
    >>> file_reader.find_between("The [deer] jumps over the fence.", '[', '[')
    'deer] jumps over the fence.'
    """
    start_idx = text.find(start)
    if start_idx == -1:
        return ""
    start_idx += len(start)
    end_idx = text.find(end, start_idx)
    if end_idx == -1:
        return text[start_idx:]
    return text[start_idx:end_idx]

def feasst_path_replace(text: str):
    """
    If the path begins with "/feasst" then replace that with the installed site-packages package
    data directory.
    """
    search = '/feasst'
    if text[:len(search)] == search:
        return str(files('feasst').joinpath('data'+text[len(search):]))
    return text

def get_mc_stage(script: str, criteria=None):
    """
    Return the current stage of the simulation based on the script.

    stage 0: Random not set
    stage 1: Configuration not set
    stage 2: Potential not set
    stage 3: ThermoParams not set
    stage 4: Criteria not set
    stage 5: Restart or Criteria is set

    >>> from feasst import file_reader
    >>> file_reader.get_mc_stage("MonteCarlo")
    0
    >>> file_reader.get_mc_stage("MonteCarlo\\nRandomMT19937")
    1
    >>> file_reader.get_mc_stage("MonteCarlo\\nRandomMT19937\\nConfiguration")
    2
    >>> file_reader.get_mc_stage("MonteCarlo\\nRandomMT19937\\nConfiguration\\nPotential")
    3
    >>> file_reader.get_mc_stage("MonteCarlo\\nRandomMT19937\\nConfiguration\\nPotential\\nThermoParams")
    4
    >>> file_reader.get_mc_stage("MonteCarlo\\nRandomMT19937\\nConfiguration\\nPotential\\nThermoParams\\nMetropolis")
    5
    >>> file_reader.get_mc_stage("Restart")
    5
    """
    stage = 0
    if 'Random' in script:
        stage = 1
    if 'Configuration' in script:
        stage = 2
    if 'Potential' in script:
        stage = 3
    if 'ThermoParams' in script:
        stage = 4
    if not criteria:
        criteria = ['Metropolis', 'AlwaysReject', 'FlatHistogram', 'MayerSampling']
    if 'Restart' in script or any(crit in script for crit in criteria):
        stage = 5
    return stage

def get_configs(script: str):
    """
    Return the lines that being with Configuration

    >>> from feasst import file_reader
    >>> file_reader.get_configs("MonteCarlo\\nConfiguration\\nConfiguration")
    ['Configuration', 'Configuration']
    >>> file_reader.get_configs(r'''MonteCarlo
    ... Configuration name=liquid
    ... Configuration''')
    ['Configuration name=liquid', 'Configuration']
    >>> file_reader.get_configs(r'''MonteCarlo
    ... Configuration
    ... Configuration name=vapor side_length=1,1
    ... Potential Model=LennardJones''')
    ['Configuration', 'Configuration name=vapor side_length=1,1']
    """
    lines = []
    for line in script.split('\n'):
        search = 'Configuration'
        if line[:len(search)] == search:
            lines.append(line)
    return lines

def get_config_names(script: str):
    """
    Return a dictionary of the names of configurations with their descriptions.

    >>> from feasst import file_reader
    >>> file_reader.get_config_names("MonteCarlo\\nConfiguration\\nConfiguration")
    {'0': 'Configuration', '1': 'Configuration'}
    >>> file_reader.get_config_names(r'''MonteCarlo
    ... Configuration name=liquid
    ... Configuration''')
    {'liquid': 'Configuration name=liquid', '1': 'Configuration'}
    >>> file_reader.get_config_names(r'''MonteCarlo
    ... Configuration
    ... Configuration name=vapor side_length=1,1
    ... Potential Model=LennardJones''')
    {'0': 'Configuration', 'vapor': 'Configuration name=vapor side_length=1,1'}
    """
    configs = {}
    for line in get_configs(script):
        name = find_between(line, start=' name=', end=' ')
        if name == "":
            name = str(len(configs))
        configs[name] = line
    return configs

def get_groups(script: str):
    """
    Return a dictionary of the names of groups with their descriptions.

    >>> from feasst import file_reader
    >>> file_reader.get_groups(r'''MonteCarlo
    ... Configuration group=fluid fluid_particle_type=0 fluid_site_type=0''')
    {'fluid': 'Configuration:0 particle_type=0 site_type=0'}
    >>> file_reader.get_groups(r'''MonteCarlo
    ... Configuration
    ... Configuration group=fluid,mof fluid_particle_type=0 mof_particle_type=1
    ... Configuration group=water water_particle_type=2 name=vapor''')
    {'fluid': 'Configuration:1 particle_type=0', 'mof': 'Configuration:1 particle_type=1', 'water': 'Configuration:vapor particle_type=2'}
    """
    config_names = list(get_config_names(script))
    groups = {}
    for iline, line in enumerate(get_configs(script)):
        grps = find_between(line, start=' group=', end=' ')
        for grp in grps.split(','):
            if grp != "":
                groups[grp] = 'Configuration:'+config_names[iline]
                ln = copy.deepcopy(line)
                #print('ln', ln)
                #print(ln.find(grp+'_'))
                for i in range(int(1e3)):
                    idx = ln.find(grp+'_')
                    #print('idx', idx)
                    if idx == -1:
                        break
                    ln = ln[idx+len(grp)+1:]
                    #print('ln', ln)
                    end = ln.find(' ')
                    if end == -1:
                        groups[grp] += ' '+ln
                    else:
                        groups[grp] += ' '+ln[:end]
    return groups

def get_particles(script: str):
    """
    Return a dictionary of the names of particles with their base file name and full path.

    >>> from feasst import file_reader
    >>> file_reader.get_particles(r'''MonteCarlo
    ... Configuration cubic_box_length=8 particle_type=lj_new.txt,atom_new.txt''')
    ({'0': 'lj_new.txt', '1': 'atom_new.txt'}, {'0': 'lj_new.txt', '1': 'atom_new.txt'})
    >>> file_reader.get_particles(r'''MonteCarlo
    ... Configuration cubic_box_length=8 particle_type=fluid:/feasst/particle/lj_new.txt''')
    ({'fluid': 'lj_new.txt'}, {'fluid': '/feasst/particle/lj_new.txt'})
    """
    particles = {}
    partdirs = {}
    for line in get_configs(script):
        arg = find_between(line, start=' particle_type=', end=' ')
        for pt in arg.split(','):
            pair=pt.split(':')
            if len(pair) == 1:
                name = str(len(particles))
                file = pair[0]
            else:
                name = pair[0]
                file = pair[1]
            idx = file.rfind('/') # remove directories
            if idx == -1:
                particles[name] = file
            else:
                particles[name] = file[idx+1:]
            partdirs[name] = file
    return (particles, partdirs)

def get_lines_in_section(particle_file: str, section_title: str):
    """
    Return the lines in a FEASST particle file with the given section_title.
    For any particle_file beginning with /feasst, replace with site-packages data directory.

    >>> from feasst import file_reader
    >>> file_reader.get_lines_in_section("/feasst/particle/lj_new.txt", 'Site Types')
    ['LJ sigma=1.0 epsilon=1.0 cutoff=3.0']
    >>> file_reader.get_lines_in_section("../../particle/spce_new.txt", 'Site Types')
    ['O cutoff=10 cutoff_outer=1 charge=-0.8476 sigma=3.16555789 epsilon=0.650169581', 'H cutoff=10 cutoff_outer=1 charge=0.4238 sigma=0 epsilon=0']
    """
    section_lines = []
    particle_file = feasst_path_replace(particle_file)
    with open(particle_file, 'r') as file1:
        lines = file1.read().splitlines()
    section_found = section_read = False
    for line in lines:
        #print('line', line, 'section_found', section_found)
        if line == "":
            if section_read:
                break
        elif line == section_title:
            section_found = True
        elif section_found:
            section_lines.append(line)
            section_read = True
    return section_lines

def get_sites_in_particle(particle_file: str):
    """
    Read the site type names from a particle file

    >>> from feasst import file_reader
    >>> file_reader.get_sites_in_particle("../../particle/lj_new.txt")
    {'site_types': ['LJ'], 'params': ['sigma', 'epsilon', 'cutoff'], 'sites': ['LJ1'], 'bonds': [], 'angles': [], 'dihedrals': []}
    >>> file_reader.get_sites_in_particle("../../particle/spce_new.txt")
    {'site_types': ['O', 'H'], 'params': ['cutoff', 'cutoff_outer', 'charge', 'sigma', 'epsilon'], 'sites': ['O1', 'H1', 'H2'], 'bonds': ['O1H1', 'O1H2'], 'angles': ['H1O1H2'], 'dihedrals': []}
    """
    sites = {'site_types': [], 'params': [], 'sites': [], 'bonds':[], 'angles': [], 'dihedrals': []}
    site_props_found = site_props_read = sites_found = sites_read = False
    for line in get_lines_in_section(particle_file, 'Site Types'):
        idx=line.find(' ')
        sites['site_types'].append(line[:idx])
        # read params in the "[site] [param]=[value]" syntax
        for pair in line[idx+1:].split(' '):
            param = pair.split('=')
            if param[0] not in sites['params']:
                sites['params'].append(param[0])
    for line in get_lines_in_section(particle_file, 'Sites'):
        idx=line.find(' ')
        sites['sites'].append(line[:idx])
    for line in get_lines_in_section(particle_file, 'Bonds'):
        idx=line.find(' ')
        sites['bonds'].append(line[:idx])
    for line in get_lines_in_section(particle_file, 'Angles'):
        idx=line.find(' ')
        sites['angles'].append(line[:idx])
    for line in get_lines_in_section(particle_file, 'Dihedrals'):
        idx=line.find(' ')
        sites['dihedrals'].append(line[:idx])
    return sites

def get_sites(script: str):
    """
    Return the sites from all particle_files in the script.
    If duplicates are present, then add a '-' to the name (except for params)
    Also, second argument in tuple return is in the same structure as the first, but includes the
    path to the files which contained those sites, etc, to include in the description when selecting.

    >>> from feasst import file_reader
    >>> file_reader.get_sites(r'''MonteCarlo
    ... Configuration cubic_box_length=8 particle_type=lj:../../particle/lj_new.txt,spce:../../particle/spce_new.txt''')[0]
    {'site_types': ['LJ', 'O', 'H'], 'params': ['sigma', 'epsilon', 'cutoff', 'cutoff_outer', 'charge'], 'sites': ['LJ1', 'O1', 'H1', 'H2'], 'bonds': ['O1H1', 'O1H2'], 'angles': ['H1O1H2'], 'dihedrals': []}
    >>> file_reader.get_sites(r'''MonteCarlo
    ... Configuration cubic_box_length=8 particle_type=lj1:../../particle/lj_new.txt,lj2:../../particle/lj_new.txt,lj3:../../particle/lj_new.txt''')[0]
    {'site_types': ['LJ', 'LJ-', 'LJ--'], 'params': ['sigma', 'epsilon', 'cutoff'], 'sites': ['LJ1', 'LJ1-', 'LJ1--'], 'bonds': [], 'angles': [], 'dihedrals': []}
    """
    particles = get_particles(script)
    sites = {'site_types': [], 'params': [], 'sites': [], 'bonds':[], 'angles': [], 'dihedrals': []}
    site_des = copy.deepcopy(sites)
    for part in particles[0]:
        pfile = particles[1][part]
        tsites = get_sites_in_particle(pfile)
        for parm in sites:
            dups = list(set(tsites[parm]) & set(sites[parm]))
            attempt = 0
            while len(dups) > 0:
                for dup in dups:
                    idx = tsites[parm].index(dup)
                    assert idx != -1
                    tsites[parm].pop(idx)
                    if parm != 'params':
                        tsites[parm] += [dup+'-']
                dups = list(set(tsites[parm]) & set(sites[parm]))
                attempt += 1
                assert attempt < 1e3
            sites[parm] += tsites[parm]
            if parm != 'params':
                site_des[parm] += len(tsites[parm])*[pfile]
    return (sites, site_des)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
