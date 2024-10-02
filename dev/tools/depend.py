"""
This utility automatically generates two kinds of files.
The first are all of the plugin/[module]/doc rst files for documentation.
The second is the feasst.h file used for the C++ interface.
In previous versions, this use to also generate the swig python interface file, py/feasst.i.
usage: /path/to/feasst/py/run.sh /path/to/feasst/dev/tools/depend.py -s /path/to/feasst
"""

import os
import sys
import re
import argparse

def parse():
    """ Parse command line arguments """
    parser = argparse.ArgumentParser()
    opt = parser.add_argument_group('opt arguments')
    opt.add_argument("--update_doc", "-u", type=int, default=1,
                     help="update documentation rst?")
    opt.add_argument("--check_headers", "-c", type=str, default='OFF',
                     help="if ON, make a cpp that includes each header to check self-sufficiency")
    required = parser.add_argument_group('required arguments')
    required.add_argument("--source_dir", "-s", help="/path/to/feasst", type=str, required=True)
    args = parser.parse_args()
    print('args', args)
    return args

def is_header(file_name):
    """ Return true if file_name is a header file """
    if file_name.endswith(".h"):
        return True
    return False

def included(file_name):
    """ Extract the included header files from a header """
    includes = list()
    if is_header(file_name):
        with open(file_name, 'r') as fle:
            lines = fle.readlines()
            for line in lines:
                if line[0:10] == '#include "':
                    includes.append(line.split('"')[1::2][0])
    return includes

def dependency(args, include_plugin, verbose = False):
    path = args.source_dir + '/plugin/'
    """ Return a list of all header files and their included files """
    depends = list()
    headers = list()
    for dir_, _, files in os.walk(path):
        if [d for d in include_plugin if d in dir_]:
            for file_name in files:
                rel_dir = os.path.relpath(dir_, path)
                rel_file = os.path.join(rel_dir, file_name)
                if '/include/' in rel_file and [d for d in include_plugin if d in rel_file]:
                    if is_header(rel_file):
                        headers.append(rel_file)
                        depends.append([rel_file, included(path+rel_file)])
                        if verbose:
                            print('***************************')
                            print(rel_file, included(path+rel_file))

    # optionally print self-sufficient check of header files
    if 'check_headers' in args:
        if args.check_headers == 'ON':
            for header in headers:
                print(header)
                #print(header[:-2])
                cpp = args.source_dir + '/plugin/' + header[:-2] + '_tmphcheck.cpp'
                cpp = cpp.replace('/include/', '/src/')
                print(cpp)
                with open(cpp, 'w') as file1:
                    file1.write('#include "'+header+'"')
    return depends

def bubble_sort(depends, verbose):
    """ bubble sort by moving headers up if they include a file not below them. """
    iteration = 0
    bubbling = True
    while bubbling:
        iteration += 1
        assert iteration < 1e4  # something is wrong with the headers
        prev = list()
        bubbling = False
        if verbose:
            print('************************')
            print(depends)
        for index, dep in enumerate(depends):
            prev.append(dep[0])
            buble = 0
            last = False
            if index == len(depends) - 1:
                last = True
            for idep in dep[1]:
                if idep not in prev:
                    buble = 1
                    if verbose:
                        print('idep', idep)
            if buble == 1 and not last:
                depends[index], depends[index+1] = depends[index+1], depends[index]
                prev[-1] = depends[index][0]
                bubbling = True
    return depends

def read_class(file_name):
    """ read the classes from a file """
    cls = list()
    with open(file_name, 'r') as fle:
        lines = fle.readlines()
        for line in lines:
            if re.search(r'^class ', line) and not re.search(';$', line):
                cls.append(line.split(' ')[1])
    return cls

def write_cpp(plugin_dir, deps, include_plugin):
    """" Write the C++ interface feasst.h """
    with open(plugin_dir+'feasst/include/feasst.h', 'w') as fsth:
        for dep in deps:
            fsth.write('#include \"' + dep[0] + '\"\n')
        # create an object in each plugin to force initialization of deserialize_map
        select_classes = list()
        plugin_classes = [
            ['beta_expanded', "ComputeBeta"],
            ['chain', "TrialGrow"],
            ['chain', "TrialParticlePivot"],
            #['cluster', "TrialRigidCluster"],
            #['configuration', ""],
            ['confinement', "ModelHardShape"],
            ['shape', "Sphere"],
            ['egce', "AEqualB"],
            ['charge', "Ewald"],
            ['example', "ModelExample"],
            #['flat_histogram', ""],
            #['math', ""],
            ['mayer', "MayerSampling"],
            ['models', "LennardJonesAlpha"],
            #['monte_carlo', ""],
            ['morph', "ComputeMorph"],
            #['opt_lj', ""],
            ['patch', "VisitModelInnerPatch"],
            ['aniso', "VisitModelInnerTable"],
            ['gibbs', "ComputeGibbsParticleTransfer"],
            ['model_expanded', "MacrostateModel"],
            ['server', "Listen"],
            ['mpi', "MPIPlaceHolder"],
            #['prefetch', ""],
            #['', ""],
            ['steppers', "Tune"],
            ['fftw', "ScatteringFFTW"],
            ['netcdf', "FileNETCDF"],
        ]
        for plugin_class in plugin_classes:
            if plugin_class[0] in include_plugin:
                select_classes.append(plugin_class[1])
        for icl in select_classes:
            fsth.write("std::shared_ptr<feasst::" + icl + "> __feasst__" + icl + " = std::make_shared<feasst::" + icl + ">();\n")

def write_docs(plugin_dir, deps, include_plugin, args, verbose):
    """ Write the docs """
    if args.update_doc == 0:
        sys.exit()
    doc = ''
    for mod in next(os.walk(plugin_dir))[1]:
        if [d for d in include_plugin if d in mod]:
            if verbose:
                print('mod', mod, 'cd', args.source_dir+'/plugin/')
            with open(plugin_dir+mod+'/doc/toc.rst', 'w') as toc:
                toc.write('\n.. toctree::\n\n')
                for dep in deps:
                    header = dep[0]
                    cls = read_class(plugin_dir+header)
                    if verbose:
                        print('cls', cls, 'header', header)
                    if re.search(mod+'/include/', header):
                        if cls:
                            doc = mod + '/doc/' + cls[0]
                            #doc = mod + '/doc/' + cls[0] + '.rst'
                            toc.write('   ' + cls[0] + '\n')
                        else:
                            funcfile = re.sub(mod+r'/include/', '', re.sub(r'.h$', '', header))
                            if verbose:
                                print('functfile', funcfile, 'mod', mod)
                                print('funcfile', funcfile)
                            toc.write('   ' + funcfile + '\n')
                            if 'include' in funcfile:
                                doc = re.sub(r'include', r'doc', funcfile)
    #                            if 'rst' not in funcfile:
    #                              doc = doc + '.rst'
                            else:
                                #funcfile[:len(mod)] == mod:
                                doc = mod + '/doc/' + funcfile
                                #doc = mod + '/doc/' + funcfile + '.rst'
                                funcfile = mod + '/include/' + funcfile
                            if verbose:
                                print('doc', doc)
                        for arguments in [True, False]:
                            extra = ''
                            if arguments:
                                extra = '_arguments'
                            with open(plugin_dir+doc+extra+'.rst', 'w') as fle:
                                if arguments:
                                    extra = ':membergroups: Arguments'
                                if verbose:
                                    print('doc', doc)
                                if cls:
                                    fle.write(cls[0]+'\n')
                                else:
                                    fle.write(funcfile + '\n')
                                fle.write('=====================================================\n')
                                if cls:
                                    for clas in cls:
                                        fle.write('\n')
                                        fle.write('.. doxygenclass:: feasst::' + clas + '\n')
                                        fle.write('   :project: FEASST\n')
                                        fle.write('   :members:\n')
                                        fle.write('   '+extra+'\n')
                                else:
                                    fle.write('\n')
                                    fle.write('.. doxygenfile:: ' + funcfile + '.h\n   :project: FEASST\n')

def read_plugins(filename):
    """ Read plugins.txt written by CMakeLists.txt """
    with open(filename, "r") as plugin_file:
        include_plugin = plugin_file.readlines()[0].split(' ')
    return include_plugin

if __name__ == '__main__':
    ARGS = parse()
    PLUGIN_DIR = ARGS.source_dir+'/plugin/'
    VERBOSE = False
    #VERBOSE = True
    INCLUDE_PLUGIN = read_plugins("plugins.txt")
    print(INCLUDE_PLUGIN)
    DEPS = dependency(ARGS, INCLUDE_PLUGIN, verbose=VERBOSE)
    DEPS = bubble_sort(DEPS, VERBOSE)
    write_cpp(PLUGIN_DIR, DEPS, INCLUDE_PLUGIN)
    write_docs(PLUGIN_DIR, DEPS, INCLUDE_PLUGIN, ARGS, VERBOSE)
