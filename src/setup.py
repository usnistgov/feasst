#! /usr/bin/env python
"""setup.py file for SWIG python interface"""
import re, os, sys
from glob import glob
from distutils.core import setup, Extension
import multiprocessing

nproc=multiprocessing.cpu_count()
print("nproc", nproc)

# read in source directory and version from command line 
#  as the 4th and 5th options, respectively
print sys.argv
source=sys.argv[3]
version=sys.argv[4]
sys.argv.remove(source)
sys.argv.remove(version)
print(source, version)

#os.environ["CC"] = "mpicxx"
#os.environ["CXX"] = "mpicxx"
#os.environ["CC"] = "icpc"
#os.environ["CXX"] = "icpc"
os.environ["CC"] = "ccache g++"   #ccache speeds up compilation
os.environ["CXX"] = "g++"
home = os.getenv("HOME")

# set compiler preprocessor macros
macs=["-DCPLUSPLUS","-DXDRFILE_H_","-DOMP_H_","-DFFTW_","-DJSON_",
      "-DFEASST_SRC_="+source, "-DVERSION="+version]

# simple utility function to search and replace list of strings add _wrap.cxx
cargs=["-fPIC","-O3","-std=c++11","-I"+home+"/include/xdrfile",
       "-L"+home+"/lib","-fopenmp","-lxdrfile","-lfftw3",
       "-I"+home+"/software/fftw-3.3.4/build/include",
       "-L"+home+"/software/fftw-3.3.4/build/lib"]+macs
largs=cargs

# monkey-patch for parallel compilation
def parallelCCompile(self, sources, output_dir=None, macros=None, include_dirs=None, debug=0, extra_preargs=None, extra_postargs=None, depends=None):
    # those lines are copied from distutils.ccompiler.CCompiler directly
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)
    # parallel code
    N=nproc # number of parallel compilations
    import multiprocessing.pool
    def _single_compile(obj):
        try: src, ext = build[obj]
        except KeyError: return
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)
    # convert to list, imap is evaluated on-demand
    list(multiprocessing.pool.ThreadPool(N).imap(_single_compile,objects))
    return objects
import distutils.ccompiler
if (nproc > 1): distutils.ccompiler.CCompiler.compile=parallelCCompile

# obtain source files by searching for *.i files in src 
#  and replacing with ".cc" and "_wrap.cxx" 
interfaces = glob("*.i")
src = list(interfaces)
#src = list()
for i in interfaces:
  #j = re.sub(r'\.i', "_wrap.cpp", i)
  #src.append(str(j))
  if (i != "feasst.i"):
    j = re.sub(r'\.i', ".cc", i)
    src.append(str(j))

feasst_module = Extension('_feasst',sources=src, extra_compile_args=cargs, 
  language="c++", extra_link_args=largs, swig_opts = ["-c++","-threads"]+macs)

setup (name = 'feasst',
       version = '0.1',
       author = "Harold Wickes Hatch",
       description = """SWIG python interface to C++ source code""",
       ext_modules = [feasst_module],
       py_modules = ["feasst"],
       )

