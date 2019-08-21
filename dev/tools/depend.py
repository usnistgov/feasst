'''
  This utility automatically generates the swig interface file.
  It also generates some of the rst files needed for documentation.
'''
usage='/path/to/feasst/py/run.sh /path/to/feasst/dev/tools/depend.py -s /path/to/feasst'

import os
import argparse
parser = argparse.ArgumentParser()
required = parser.add_argument_group('required arguments')
required.add_argument("--source_dir", "-s", help="/path/to/feasst", type=str, required=True)
args = parser.parse_args()
import sys
import pyfeasst

external_libs=['xdrfile.h', 'xdrfile_xtc.h']
#include_plugin=['utils','math','configuration','system','steppers','monte_carlo','flat_histogram']
include_plugin=['utils','math','configuration','system','steppers','monte_carlo','mayer','models','patch','example','confinement','chain','flat_histogram','ewald','growth_expanded']

def is_header(file_name):
  if file_name.endswith(".h"):
    return True
  return False

# extract the included header files from a header
def included(file_name):
  includes = list()
  if is_header(file_name):
    with open(file_name, 'r') as fle:
      lines = fle.readlines()
      for line in lines:
        if line[0:10] == '#include "':
          includes.append(line.split('"')[1::2][0])
  return includes

# return a list of all header files and their included files
def dependency(path):
  depends = list()
  headers = list()
  for dir_, _, files in os.walk(path):
    if [d for d in include_plugin if d in dir_]:
      for fileName in files:
        relDir = os.path.relpath(dir_, path)
        relFile = os.path.join(relDir, fileName)
        if '/include/' in relFile and [d for d in include_plugin if d in relFile]:
          if is_header(relFile):
            headers.append(relFile)
            depends.append([relFile, included(relFile)])
            #print('***************************')
            #print(relFile, included(relFile))
  # check if the includes are clean
  for dep in depends:
    for dep1 in dep[1]:
      if dep1 not in (headers + external_libs):
        raise Exception(dep1, 'is included by', dep[0], 'but has wrong directory structure')
  return depends

# bubble sort by moving headers up if they include a file not below them.
def bubble_sort(depends):
  iteration = 0
  bubbling = True
  while bubbling:
    iteration += 1
    assert(iteration < 1e4) # something is wrong with the headers
    prev = list()
    bubbling = False
    #print('************************')
    #print(depends)
    for index, dep in enumerate(depends):
      prev.append(dep[0])
      buble = 0
      last = False
      if index == len(depends) - 1:
        last = True
      for idep in dep[1]:
        if idep not in prev:
          buble = 1
          #print('idep',idep)
      if buble == 1 and not last:
        depends[index], depends[index+1] = depends[index+1], depends[index]
        prev[-1] = depends[index][0]
        bubbling = True
  return depends

# read the classes from a file
import re
def read_class(file_name):
  cls = list()
  with open(file_name, 'r') as fle:
    lines = fle.readlines()
    for line in lines:
      if re.search(r'^class ', line) and not re.search(';$', line):
        cls.append(line.split(' ')[1])
  return cls

# obtain the headers sorted by dependency
# write the swig interface file
with pyfeasst.cd(args.source_dir+'/plugin/'):
  deps=bubble_sort(dependency('./'))
  classes=list()
  # note this method for reading classes is copied below for doc
  for dep in deps:
    header = dep[0]
    cls = read_class(header)
    classes.append(cls)

# obtain the headers sorted by dependency
# write the swig interface file
with pyfeasst.cd(args.source_dir+'/plugin/'):
  with open('../py/feasst.i', 'w') as swig_file:
    swig_file.write(
"/* This is an interface file for swig.\n\
   This file is created by dev/tools/depend.py . Modifications to this\n\
   file will likely be overwritten in the future. Thus, edit depend.py\n\
   instead.\n\
\n\
   usage: "+usage+"\n\
 */\n\
\n\
%module feasst\n\
\n\
%{\n\
")
    for dep in deps:
      swig_file.write('#include "' + dep[0] + '"\n')
    swig_file.write(
"using namespace feasst;\n%}\n\
\n\
%include \"std_string.i\"\n\
%include \"std_vector.i\"\n\
%include \"std_shared_ptr.i\"\n\
%include \"std_iostream.i\"\n\
%template(IntVector) std::vector<int>;\n\
%template(IntIntVector) std::vector<std::vector<int> >;\n\
%template(DoubleVector) std::vector<double>;\n\
%template(AnalyzeVector) std::vector<std::shared_ptr<Analyze> >;\n\
using namespace std;\n\
%pythonnondynamic;\n\
%include \"std_map.i\"\n\
%template(args) std::map<std::string, std::string>;\n\
\n\
")
    for cls in classes:
      for icl in cls:
        swig_file.write("%shared_ptr(feasst::" + icl + ");\n")
    for dep in deps:
      swig_file.write('%include ' + dep[0] + '\n')

# write the docs
with pyfeasst.cd(args.source_dir+'/plugin/'):
  for mod in next(os.walk('.'))[1]:
    if [d for d in include_plugin if d in mod]:
      print('mod', mod, 'cd', args.source_dir+'/plugin/')
      with open(mod+'/doc/toc.rst', 'w') as toc:
        toc.write('API\n************\n\n.. toctree::\n\n')
        for dep in deps:
          header = dep[0]
          cls = read_class(header)
          print('cls', cls, 'header', header)
          if re.search(mod, header):
            if cls:
              if 'utils' in doc:
                print("hi", doc)
              doc = mod + '/doc/' + cls[0] + '.rst'
              toc.write('   ' + cls[0] + '\n')
            else:
              funcfile = re.sub(mod+r'/include/', '', re.sub(r'.h$', '', header))
              print('functfile', funcfile, 'mod', mod)
              print('funcfile', funcfile)
              toc.write('   ' + funcfile + '\n')
              if 'include' in funcfile:
                doc = re.sub(r'include', r'doc', funcfile)
                if 'rst' not in funcfile:
                  doc = doc + '.rst'
              else:
                #funcfile[:len(mod)] == mod:
                doc = mod + '/doc/' + funcfile + '.rst'
                funcfile = mod + '/include/' + funcfile
              print('doc', doc)
            with open(doc, 'w') as fle:
              print('doc', doc)
              if cls:
                fle.write(cls[0]+'\n')
              else:
                fle.write(funcfile + '\n')
              fle.write('=====================================================\n')
              if cls:
                for cl in cls:
                  fle.write('\n')
                  fle.write('.. doxygenclass:: feasst::' + cl + '\n')
                  fle.write('   :project: FEASST\n')
                  fle.write('   :members:\n')
              else:
                  fle.write('\n')
                  fle.write('.. doxygenfile:: ' + funcfile + '.h\n   :project: FEASST\n')
