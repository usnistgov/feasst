# FEASST
>The Free Energy and Advanced Sampling Simulation Toolkit (FEASST) is a free, 
open-source, modular program to conduct molecular and particle-based 
simulations with flat-histogram Monte Carlo and molecular dynamics methods.

### Features
- Wang-Landau, Transition-Matrix and/or Metropolis Monte Carlo
- Canonical, grand canonical and expanded ensembles
- Interface as a python module or C++ class
- Easy use of MPI and checkpointing
- Many advanced sampling Monte Carlo moves
- Supported energy functions include anisotropic particles and the Ewald
summation
- Robust unit testing

### Advanced Monte Carlo moves
- Configurational bias insertions, deletions and regrowth with multiple first
bead insertion
- Aggregation volume bias (AVB) insertions, deletions and the AVB2  and AVB3 
algorithms
- Geometric cluster algorithm
- Rigid cluster moves
- Parallel configuration swaps
- Constant pressure and floppy box

### Intermolecular interactions
- Charged interactions with the Ewald summation
- Lennard Jones with Yukawa, LRC, force shift, Gaussians
- Patchy particles
- Hard spheres/square well
- Tabular potentials for anisotropic particles

# Contents
* [Usage](#usage)
* [Prerequisites](#prerequisites)
* [Installation](#installation)
* [External libraries](#external-libraries)
* [Test case 1. Lennard-Jones](#test-case-1-lennard-jones)
* [Analysis of configurations for WL-TMMC simulations](#analysis-of-configurations-for-wl-tmmc-simulations)
* [Initializing a simulation from an XYZ file](#initializing-a-simulation-from-an-xyz-file)
* [Example of custom analysis in input script](#example-of-custom-analysis-in-input-script)
* [Example of restarting a simulation](#example-of-restarting-a-simulation)
* [Example of an isotropic tabular potential](#example-of-an-isotropic-tabular-potential)
* [Example of adding or modifying a pair potential](#example-of-adding-or-modifying-a-pair-potential)
* [Example of adding or modifying a Monte Carlo trial](#example-of-adding-or-modifying-a-monte-carlo-trial)
* [Example of adding or modifying an analysis code](#example-of-adding-or-modifying-an-analysis-code)
* [Example of adding or modifying a random number generator](#example-of-adding-or-modifying-a-random-number-generator)
* [Contact](#contact)

# Usage

The following example may be found in [example/lj](example/lj).

In C++, a simple NVT Lennard Jones (LJ) simulation is performed as follows:
```c++
#include "pair_lj.h"
#include "mc.h"
#include "ui_abbreviated.h"

int main() {
  feasst::Space space(3, 0);
  for (int dim = 0; dim < space.dimen(); ++dim) {
    space.lset(8, dim); // 8 box length
  }
  space.addMolInit("../../forcefield/data.lj");
  feasst::PairLJ pair(&space, 3);   // potential truncation at 3
  pair.initEnergy();
  feasst::CriteriaMetropolis criteria(1.2, 1.);  // 1/kT = 1.2
  feasst::MC mc(&space, &pair, &criteria);
  feasst::transformTrial(&mc, "translate", 0.1);
  mc.nMolSeek(50, "../../forcefield/data.lj");  // add 50 particles
  mc.initLog("log", 1e4);
  mc.initMovie("movie", 1e4);
  mc.runNumTrials(1e6);
}
```

This simulation is compiled and run by a bash script [example/lj/run.sh](example/lj/run.sh):
```bash
prog=nvt
$FEASST_INSTALL_DIR_/tools/compile.sh $prog
./$prog
```
Note that `FEASST_INSTALL_DIR_` should be specified in your `~/.bash_profile` as `export FEASST_INSTALL_DIR_="$HOME/path/to/feasst/"`.

In python, the same simulation may be written as in the file [example/lj/lj.py](example/lj/lj.py)
```python
#! /usr/bin/env python
import os, sys
feasstdir = os.getenv("FEASST_INSTALL_DIR_") + "/build"
if (not os.path.isfile(feasstdir+"/_feasst.so")):
  feasstdir = os.getenv("FEASST_INSTALL_DIR_") + "/src"
sys.path.append(feasstdir)
import feasst
space = feasst.Space(3, 0)
for dim in range(space.dimen()): space.lset(8, dim) # 8 box length
space.addMolInit("../../forcefield/data.lj")
pair = feasst.PairLJ(space, 3)    # potential truncation at 3
pair.initEnergy()
criteria = feasst.CriteriaMetropolis(1.2, 1.);  # 1/kT = 1.2
mc = feasst.MC(space, pair, criteria)
maxMoveParam = 0.1
feasst.transformTrial(mc, "translate", maxMoveParam)
mc.nMolSeek(50, "../../forcefield/data.lj")   # add 50 particles
mc.initLog("log", int(1e4))
mc.initMovie("movie", int(1e4))
mc.runNumTrials(int(1e6))
```

And for those who prefer a text-based input file, see [example/lj/text/input.txt](example/lj/text/input.txt):

```
space 3   # simulations in 3D
boxl 8    # 8 box length
addMolInit ../../forcefield/data.lj
pair lj 3   # rCut of 3
initEnergy
criteria metropolis 1.2 1.  # 1/kT = 1.2
mc
trial translate 0.1   # 0.1 max move parameter
nMolSeek 50 ../../forcefield/data.lj
log log 10000
movie movie 10000
run 1000000
```

This text-based simulation is then run as `/path/to/feasst/[build/bin,src]/ui_text -f input.txt`.
Commands are interpreted via [src/ui_text.cc](src/ui_text.cc), which is cumbersome to maintain.
It is highly recommended to use the C++ or python interface instead.

# Prerequisites

FEASST is designed for a LINUX or MAC OS X platform with the following minimum version software.

- make >= 3.81
- compiler with c++0x support (e.g., g++ >= 4.7)

### Optional tools:
- CMake >= 2.8.12.2
- SWIG >= 1.3.40 (python interface)
- anaconda >= 1.9.1 (python >= 2.7)
- xdrfile 1.1b (compressed xtc trajectories)
- gtest >= 1.7.0 (C++ unittests)
- valgrind (C++ memory testing for development)
- doxygen >= 1.6.1 (C++ documentation)
- openmpi >= 1.4.5 (parallel computation)

# Installation

Installation may be performed with either CMake or a plain Makefile.
If you do not have a preference, it is recommended to attempt CMake.
Otherwise, the plain Makefile approach is available.
The example input scripts automatically check for the CMake install first,
so make sure that you remove your CMake build files if you want to use
the install from the [src/Makefile](src/Makefile) instead.

## CMake

CMake is the recommended installation method.

```bash
cp -r buildtemplate build
cd build
cmake .
make -j 8
```

To control the install, you can edit `CMakeLists.txt` in `build` as follows 
before running the `cmake .` command.

For example, to use the XDRFILE library for xtc files:
```cmake
option(USE_XDRFILE "Use xdrfile library" OFF)
```
Or
```bash
cmake -DUSE_XDRFILE=ON .
```

Or to give CMake the path to your xdrfile library:
```cmake
set(XDRFILE_DIR "/path/to/xdrfile")
```
Or
```bash
cmake -DXDRFILE_DIR=/path/to/xdrfile .
```

Or, for example, if you want to use the python interface, then use the
following CMake command instead: 
```bash
cmake -DUSE_SWIG=ON .
```

If you are changing the default build options in `CMakeLists.txt`,
make sure to start compilation with a fresh `build` directory before CMake is
invoked (e.g., completely remove the build directory and start over, after 
saving any relevant changes to `CMakeLists.txt`).

## Makefile

This is the old compilation method, but may still be used if you prefer.

### Python installation
```bash
cd src
make swig
```

### C++ installation
```bash
cd src
make cnotest
```

### Text-based installation (using C++ engine)
```bash
cd src
make ui_text
```

For any interface, modify [src/Makefile](src/Makefile) to control external libraries (below).

# External libraries

Here is how to set up external libraries you may want to use with FEASST.
To begin, some libraries require installation. And some require certain compiler flags if not using CMake.
If you do not wish to use the libraries and CMake, make sure the compiler flags are not included in [src/Makefile](src/Makefile).

### XTC 1.1b:
For writing compressed XTC trajectory files.
```bash
ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.tar.gz
tar -xf xdrfile-1.1.tar.gz; cd xdrfile-1-1b
./configure --enable-shared --prefix=$HOME/ #enable-shared for SWIG
make install
```
Associated compiler flags in [src/Makefile](src/Makefile):
```
-DXDRFILE_H_
-DCPLUSPLUS
-I/path/to/install/dir/include/xdrfile
-L/path/to/install/dir/lib
-lxdrfile
```
### Google Test 1.7.0
For testing the C++ code:

- download gtest: http://code.google.com/p/googletest/downloads/detail?name=gtest-1.7.0.zip
- unzip and correct path in `GTEST_DIR` in [src/Makefile](src/Makefile). No need to compile gtest
- to compile the FEASST unittests, use the command `make c` in [src](src)

### OpenMPI with Intel compilers
```bash
tar -xf openmpi*gz; cd openmpi*; mkdir build; cd build
../configure --prefix=`pwd`/.. CC=icc CXX=icpc $intel compilers
make
make install
```
Associated compiler flags in [src/Makefile](src/Makefile):
```
-DMPI_H_
```

### OpenMP
Associated compiler flags in [src/Makefile](src/Makefile):
```
-DOMP_H_
-fopenmp
```

### FFTW 3.3.4
This library is used for computing the scattering of anisotropic shapes.
```bash
# download fftw-3.3.4, uncompress, move to main directory
./configure --prefix=/path/to/install/dir --enable-shared --with-pic
make
make install
```
Associated compiler flags in [src/Makefile](src/Makefile):
```
-DFFTW_
-I/path/to/install/dir/fftw-3.3.4/build/include
-L/path/to/install/dir/fftw-3.3.4/build/lib
-lfftw3
```

### VMD 1.9.2
VMD is great for visualizing and analyzing trajectories.
```bash
# download vmd
tar -xf vmd-1.9.2.bin.LINUXAMD64-RHEL5.opengl.tar.gz
cd vmd-1.9.2
# edit the configure file to change install location
./configure LINUXAMD64
make install -j 8
# add VMD to your path
PATH=$PATH:/path/to/install/dir/vmd-1.9.2/bin/
```

### SWIG 2.0.12
Required for python installation.
```bash
cd swig-2.0.12; ./configure --prefix=/path/to/install/dir; make
```

### CMake 2.8.12.2
Download from https://cmake.org/files/v2.8/
`tar -xf cmake-2.8.12-rc2-Linux-i386.tar.gz`

### HDF5 1.8.18
```bash
sudo ./configure --prefix=/usr/local/hdf5 --enable-cxx
make; make check; make install; make check-install
```
Associated compiler flags in [src/Makefile](src/Makefile):
```
-DHDF5_
-I/path/to/install/dir/include
-L/path/to/install/dir/lib
-lhdf5 -lhdf5_cpp
```

### GSL 2.3
For spline interpolation.
```bash
./configure --prefix=/path/to/install/dir; make; make install
```
Associated compiler flags in [src/Makefile](src/Makefile):
```
-DGSL_
-I/path/to/install/dir/gsl-2.3/include
-L/path/to/install/dir/gsl-2.3/lib
-lgsl -lgslcblas -lm
```
# Test case 1. Lennard Jones

The first test case is to simulate the LJ potential and reproduce the
results found in the NIST SRSW for LJ, EOS-TMMC: 

http://www.nist.gov/mml/csd/informatics_research/srsw.cfm

http://mmlapps.nist.gov/srs/LJ_PURE/eostmmc.htm

Move to the test directory
`cd ../testcase/lj/srsw/eostmmc/1.5`

The file muvttmmclj.cc is the C++ interface for FEASST.
The file run.sh is an example script to compile and run muvttmmclj.
The file data.lj is a LAMMPS style description of an LJ particle.
The files `lj.msdb.*` contain the results from the SRSW.

Run the simulation `./run.sh`, or `./muvttmmclj.py`, and it will produce:

- `log` file printing information by step
- `colMat` file with the macrostate probability distribution (lnpi), potential energy and collection matrix
- `movie*` coordinate files viewable with VMD

To compare the results with the NIST SRSW, compare the following:

- colMat columns 1:2 with `lj.msdb.t150.*.p_macro.dat`
- colMat columns 1:3 with `lj.msdb.t150.*.energy.dat`

To run the test on multiple processors using OMP:
```bash
cd omp
./run.sh OR ./muvttmmclj.py
```

The parallel version has only 3 different lines
```c++
// mc.runNumTrials(npr);  // single processor run command
mc.initWindows(1);        // initialize parallel windows
mc.runNumSweeps(10, -1);  // run until number of sweeps reached
```


# Analysis of configurations for WL-TMMC simulations

For analyzing configurations to post-process simulation data (e.g., 
log and movie files), see the example code [tools/xyz2bin.cc](tools/xyz2bin.cc)

This code both splits the xyz files for a given order parameter and shows example of doing some 
analysis within the c code. It uses the checkpointing to pick up run 
variables. For example it knows how often you printed movies versus 
logs. It also picks up on the number of processors and can average 
over all of those processors for analysis.

It outputs files:
  `analysis*` for average x position of first molecule
  moviep[proc]b[bin].xyz

where proc is the processor number and bin is the order parameter
index as described by the acceptance criteria.

# Initializing a simulation from an XYZ file

The following code reads an xyz file format to input an initial 
configuration.

```c++
std::ifstream inFile("nameOfFile");
p.readxyz(inFile);  // "p" is the name of your pair object
```

If there are zero particles in the space class, then it automatically 
attempts to add the molecules based on the first molecule type 
described by the s.addMolInit function.

In cases with multiple components, this is not sufficient. So you will
want to make sure you initialize the appropriate number of molecules 
in the appropriate order. For example, if your xyz file lists A,B,A,B.
or A,A,A...,B,B,B... then you need to add these in the right order.

For example, something like the following:

```C++
for (int iMolA = 0; iMolA < nMolA; ++iMolA) {
  s.addMol("/name/of/data/file");
}
```

Then the same for B, assuming your xyz has all A listed, followed by 
all B.

After all of the s.addMol commands are performed and the xyz file is 
read, you will need to update the pair class as follows:

```c++
p.addPart();
p.initEnergy();
```

A simple test that the xyz file was read correctly is to print it and
compare:

```c++
p.printxyz("filename",1);
```

# Example of custom analysis in input script

In C++, you can define your own custom derived Analyze class inside the input
file.

First, you can define an analysis as follows which accumulates the potential
energy:
```c++
#include "analyze.h"
class AnalyzeMonkeyPatch : public Analyze {
 public:
  AnalyzeMonkeyPatch(Space *space, Pair *pair) : Analyze(space, pair) {}
  ~AnalyzeMonkeyPatch() {}
  Accumulator pe;
  void update() {
    pe.accumulate(pair_->peTot()/double(space_->nMol()));
  }
  void print() {
    cout << pe.average() << " +/- " << pe.blockStdev() << endl;
  }
};

```

After defining your new Analyze class, you may then add it to your MC class:
```c++
AnalyzeMonkeyPatch an(&space, &pair);
an.initFreq(1);
an.initPrintFreq(1e5);
mc.initAnalyze(&an);
```

This example is shown in the test case [testcase/lj/srsw/nvt-mc/lj.cc](testcase/lj/srsw/nvt-mc/lj.cc).

Note that while this example is in the spirit of a monkey patch, implementing
a monkey patch on the SWIG python objects requires editting the vtable.
In this case, it may be easier to add the custom analysis in the source directory.
See [Example of adding or modifying an analysis code](#example-of-adding-or-modifying-an-analysis-code).

# Example of restarting a simulation

Checkpoint files may be written periodically during a simulation, and
these may be used to restart a simulation. For example, see [test/binary/tee/table/tee_rst.cc](test/binary/tee/table/tee_rst.cc)

In this file, the checkpoint file is simply read and restarted in
two lines for single processor simulations:

```c++
// read checkpoint files
feasst::WLTMMC mc("tmp/rst");

// run simulation
mc.runNumTrials(npr);
```

Note that, if you are attempting to restart a simulation that was
terminated abruptly, it is possible that the checkpoint files were in
the process of writing during the termination. In this case, the files
themselves could be missing important details. If this is the case,
your simulation will likely crash upon restart or output potential
energy which is quite different from the previous value. To remedy
this situation, the checkpoint files keep a 'backup' file which ends
in ".bak", which you may use instead. If you wish to use the backup
files, then all files ending in '.bak' should replace the same files
without the '.bak' ending. Before replacing files, it is recommended
to first backup the entire tmp folder.

Note that multiprocessor simulations may take additional care to
restart correctly. If you wish to restart just one processor, you may
simply use one of the files with an appended p# (e.g. "tmp/rstp0").

If you wish to restart simulations that are independent, then an 
example may be found in [example/binary/tee/shortrange/tee_rst.cc](example/binary/tee/shortrange/tee_rst.cc)

In this file, the two lines are as follows:

```c++
// read restart file
feasst::WLTMMC mc("tmp/rst");

// run sweeps
mc.runNumSweepsRestart(100, "tmp/rst");
```

Restarting simulations that are coupled (e.g. by configuration swaps)
may require more initialization that is not currently described in 
this documentation.

# Example of an isotropic tabular potential

Instead of implementing your own pair potential in the code, you may simply make a text file with your potential. 
An example may be found in the following test directory: [test/binary/tee/table](test/binary/tee/table)

In this example, a binary LJ-lambda potential is simulated. In tee.cc
the potential is implemented with PairLJMulti, printed, and then used
to initialize PairTabular. It outputs the tables as files `tabi*j*`
which have headers like the following:

```
# tableType PairLJMulti
# tabDims 1
# dimorder0 0
# dimn0 5001
# dimmin0 0.94089999999999996
# dimmax0 1.1664000000000001
9542.2200121376991
9483.2587236908766
9424.6627162728782
```

Note that `dimn0` is the number of table elements.
Distances are shown as a function of the variable s=r^2, such that dimmax0=rCut^2=1.08^2 and dimmin0=rCutInner^2=sigFac^2=0.97^2.
For tabular potentials, r < rCutInner has infinite potential energy.

An example of utilizing the table potential (without generating) is provided in [test/binary/tee/table/tee_nogen.cc](test/binary/tee/table/tee_nogen.cc)

This file and the `tabi*j*` files may be used as templates to create
your own pair potentials.

# Example of adding or modifying a pair potential

To begin, one may create a new pair potential by simply copying an
existing version that is closest to the type of potential you are 
trying to implement. 

For example, copy the files [src/pair_lj_multi.h](src/pair_lj_multi.h) to 
a file for your new potential, such as `src/pair_lj_multi_opt.h`
and [src/pair_lj_multi.cc](src/pair_lj_multi.cc) to `src/pair_lj_multi_opt.cc`
Within these files, replace
all instances of `PairLJMulti` with `PairLJMultiOpt`, and the header 
macro `PAIR_LJ_MULTI_H_` to `PAIR_LJ_MULTI_OPT_H_`. The makefile will 
automatically compile any file name that begins with `src/pair_`.
At this stage, one may
run a simulation of the new pair class and verify that it is 
identical to the old one.

In addition, you may rather create a tabular potential and use 
[src/pair_tabular_1d](src/pair_tabular_1d) as described in the section above.

The pair potential may be called in 3 subroutines:

1. `initEnergy()` - loops through all pairs of particles. Designed to be
   slow, this is typically called every 1e5 MC steps to prevent energy drift during MC and test against optimizations such as neighbor lists and cell lists.

2. `multiPartEner*` - computes energy of only selected atoms. Designed
   to be fast, this routine is called often and should be optimized.
   The most basic test (called by `checkEnergy(tol, 1)`) is that this
   yields the same energy as `initEnergy()`. Because errors in the pair
   potential can be difficult to test, it is recommended to implement
   the pair potential in `multiPartEner*` independently from `initEnergy()`.

   Note that, in many cases, the `multiPartEner` subroutine calls more
   specific versions based on simulation parameters (e.g., dimensions
   or cutoff types). Some of these special cases are implemented in 
   the Pair base class, and thus require only a modification to a 
   subroutine named `multiPartEner*Inner` in the derived class.

3. `allPartEnerForce` - these subroutines are utilized for large cluster
   moves or simulation domain changes. In many cases, the potential is
   implemented in a subroutine named `allPartEnerForceInner`.

Note that in some cases, the forces may be calculated. These are 
sometimes used for "smart Monte Carlo" or molecular dynamics, but otherwise are unnecessary for Monte Carlo.

You have to add the new pair class to `makePair()` in [src/pair.cc](src/pair.cc), and include the header, for restart capability.

For the Swig python interface, you also must copy [src/pair_lj_multi.i](src/pair_lj_multi.i) into a new file, with the new header,
and also add the new pair header to [src/pair.i](src/pair.i) and [src/feasst.i](src/feasst.i).

# Example of adding or modifying a Monte Carlo trial

You can begin by copy/pasting and renaming an existing trial and renaming the class.
But remember, this can carry along a lot of garbage.
So copy the simplest trial that is closest to what you want to accomplish, 
and immediately remove the unnecessary pieces.
Rename all `TrialName` to `TrialNewName` and `TRIAL_NAME_H` to `TRIAL_NEW_NAME_H_`

For the python/SWIG interface, add a new `trial_*.i` file, 
and add this file to the list in [src/trial.i](src/trial.i) as well as [src/feasst.i](src/feasst.i).

For the user interface, you may also add this new trial to [src/ui_abbreviated.h](src/ui_abbreviated.h) or some other `ui_` file.

For restart capability, add this new trial to `makeTrial` in [src/trial.cc](src/trial.cc).

The MC class only knows the base Trial class. In some cases you may need to add
a feature to the base Trial class, but try to avoid this if possible.

# Example of adding or modifying an analysis code

Adding and modifying `analyze_` files are even simpler than the pair and trial examples above.
All you need to do is copy/paste the one of the existing `analyze_` files and update the class names, headers, etc, as described above.

Alternatively, you can define your custom analysis code in the input script instead of adding it to the FEASST code base.
For an example, you can see [src/analyze_unittest.cc](src/analyze_unittest.cc) for an example `AnalyzeMonkeyPatch`.
In the energy case, you would most likely define a new `void update` function which is run at the end of every `nFreq` steps.

# Example of adding or modifying a random number generator

As described above, except use [src/random_nr3.h](src/random_nr3.h) as a template.

# Contact

Project lead: Harold Wickes Hatch

www.nist.gov/people/harold-hatch

harold.hatch@nist.gov

For list of contributors, see [CONTRIBUTORS.txt](CONTRIBUTORS.txt)

# For more advanced documentation, see [DEV.txt](Dev.txt) for developers

![alt text](http://hhatch.com/classBase.png "FEASST class structure")
