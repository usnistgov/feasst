#! /usr/bin/env python
"""setup.py file for SWIG python interface"""
import re, os, sys
from glob import glob
from distutils.core import setup, Extension

os.environ["CC"] = "g++"   #ccache speeds up compilation
os.environ["CXX"] = "g++"
#os.environ["CC"] = "mpicxx"
#os.environ["CXX"] = "mpicxx"
#os.environ["CC"] = "icpc"
#os.environ["CXX"] = "icpc"

# use external library args from Makefile
cargs = sys.argv[3]
sys.argv.remove(cargs)

# add one additional command line option
cargs += " -fPIC"
cargs = cargs.split()

# monkey-patch for parallel compilation
import multiprocessing
nproc=multiprocessing.cpu_count()
print("nproc", nproc)

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
for i in interfaces:
  if (i != "feasst.i"):
    j = re.sub(r'\.i', ".cc", i)
    src.append(str(j))

feasst_module = Extension('_feasst',sources=src, extra_compile_args=cargs,
  language="c++", extra_link_args=cargs, swig_opts = ["-c++","-threads"])

setup(name = 'feasst',
      version = '0.1',
      author = "Harold Wickes Hatch",
      description = """SWIG python interface to C++ source code""",
      ext_modules = [feasst_module],
      py_modules = ["feasst"],
      )

