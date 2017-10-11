*************
Notes
*************

These are notes that I plan to format properly and add to the documentation
at a later time.
I likely wrote these in response to a particular question.

Log file and tuning
#####################

For example, for the log file below:

[hwh@pn101924 cpp]$ head -10 log
# attempts nMol pe/nMol translate maxMove TrialSwap Na Nb vol maxMove
10000 38 -0.412251 0.94403 0.1 0.301035 8 30 0.979839 0.001
20000 38 -0.418134 0.946039 0.105 0.296319 6 32 1 0.00105
30000 38 -0.496219 0.93927 0.11025 0.282662 3 35 0.996 0.0011025
40000 38 -0.447826 0.938576 0.115763 0.281589 5 33 0.989091 0.00115763
50000 38 -0.473118 0.923352 0.121551 0.277183 4 34 0.987755 0.00121551
60000 38 -0.545913 0.928427 0.127628 0.275281 3 35 0.971698 0.00127628
70000 38 -0.430846 0.922709 0.13401 0.274337 3 35 0.978632 0.0013401
80000 38 -0.393172 0.925612 0.14071 0.278066 1 37 0.98913 0.0014071
90000 38 -0.604083 0.91554 0.147746 0.281482 2 36 0.991803 0.00147746

the 4th column, "translate" is the name of the move type, and that column gives the move acceptance. The 5th column is 'maxMove' which is the maximum displacement size for the translation. It is "tuning" ever 1e4 attempts, so you can see that the max move parameter is increasing by ~5% because the acceptance is >0.25. Trial swap has no max move parameter to tune. The 9th column is the acceptance for the "vol" move, and the 10th column i sthe maxMove parameter for the volume.

You can change the default target acceptance percentage by setting the variable "trialPointer->targAcceptPer". This requires you to define a shared_ptr<TrialTransform> trialPointer and do the mc.initTrial(trialPointer) instead of defining the move via the "ui_abbreviated" interface (e.g., trialTransform(&mc, "translate", ..) give you no access to the derived class TrialTransform).

Documentation
################

pip install sphinx
pip install breathe
doxygen with GENERATE_XML
run sphinx-quickstart, enable autodoc
add something like the following to your sphinx index.rst::

    .. doxygenclass:: Nutshell
       :project: nutshell
       :members:

add the following to your sphinx conf.py
  extensions = [ "breathe" ]
  breathe_projects = {"FEASST":"../xml"}

pip install sphinx_rtd_theme

run sphinx: make html

For math in doxygen commends, use::

    \f$ latex code here \f$


Checklist Before Public Distribution
########################################

* Documentation. Examples folder. Website.
* Check licensing on mins.h or use an alternative minimization,
  and also use of LAMMPS code (jacobi)
* Stream-line and make python interface more accessible / documented.
* move checkBond from MC to analyze class for rigid particle simulations (or as part of checkE?)
* cpplint everything
* follow the google style guide (loosely)
* Move some analysis outside of space class
  - but what if things like trials use the analysis, such as clusters?
  - ""This trial requires MC to have an analyze class...""

Packaging
#####################

git tag -a <version number>

git archive master --prefix='feasst'`git describe master`'/' | gzip > feasst`git describe master`.tar.gz

FAQ
############

Testing potentials
####################

Usually when im starting a new system I throw in a few unittests to compare the 'analytical' solution for the pair interactions, which helps me sleep better at night since there is basically no other check for errors in implementing potentials (except for wrong results!). So for example if you look in pair_lj_multi_unittest.cc there are a bunch of tests like this

WCAanalytical

cg3analytical

LJYanalytical

which a lot of these use the linearShifts. To be absolutely sure you could follow the similar procedure of using "xAdd" to place two molecules precisely in a box and test the energies.

How to efficiently divide windows
*************************************

This is a common issue with a few approaches

1. You can use initWindows(nExp= variable to set the spacing based on some exponential distribution,e.g., nExp=3. makes windows even bigger o the low density side vs nExp=2.

http://feasst.hhatch.com/WLTMMC.html#_CPPv2N6WLTMMC11initWindowsEKdKi

void initWindows(const double nExp, const int nOverlap)

2. The preliminary sampling idea is hard to implement because its hard to estimate how long it will take to converge the high density region, because you don't really know when it will converge (if ever) until it finally does

3. One thing I did for some really high density simulations is I simply broke them up into multiple simulations. For example, I launched two jobs simultaneously, one with N=0 to 336, and another N=336 to N=550. Of course the low density ones finished 10x faster but then the high density windows were smaller

4. If you're trying to do a high throughput approach where each model and condition can be quite different but you don't want to have to hand pick conditions for each one, you can do a ridiculous number of windows (N=16/32 perhaps?) and terminate after some run time and only use the windows that managed to converge and throw out the rest. Still you would want to choose carefully the nMolMax to not have too much wasted processor time.

One thing to look out for is you want the processors running high density to have some kind of access to the configurations coming up from the low density windows to help convergence so its not just stuck in some glass. In that regard it may help to have the OMP configuration swaps on (TrialConfSwapOMP) but I would hesitate to increase the frequency of these swaps because they break detailed balance.

Getting to high density
*************************

See above, but also:

You can use the more MD-style approach of initializing a big box and squeezing it.

.. code-block:: c++

   c.pressureset(1000.);
   transformTrial(&mc, "lxmod", 0.001);
   transformTrial(&mc, "lymod", 0.001);
   transformTrial(&mc, "lzmod", 0.001);
   // transformTrial(&mc, "vol", 0.001);  // this one does an isotropic volume move instead of independent. You don't need both

To remove the box move after equilibration, I prefer to "scope" the mc class with the box move and just make another clean mc class, or do an shallow copy (cloneShallow) before calling the transform trial and then use that clone later for production.
Another option is mc.removeTrial(trial#);

Recompile with -fPIC
*********************

Issue: Compilation error gives "can not be used when making a shared object; recompile with -fPIC"

Solution: One of your external libraries (e.g., fftw or xdrfile) needs the flag "--enable-shared" during configuration. Or you can edit CMakeLists.txt to add "-fPIC" as follows:

SWIG_LINK_LIBRARIES(feasst ${PYTHON_LIBRARIES} ${EXTRA_LIBS} -fPIC)   # HWH: add -fPIC
#SWIG_LINK_LIBRARIES(feasst ${PYTHON_LIBRARIES} ${EXTRA_LIBS})        # HWH: old version

Unittest that sometimes fail
##############################

[ RUN      ] Trial.allmoves
Note in /home/hwh/feasst/src/functions.cc line 89: time(seed): 1496248513
id metropolis
id wltmmc
lj metropolis
lj wltmmc
onePatch metropolis
onePatch wltmmc
/home/hwh/feasst/src/trial_unittest.cc:297: Failure
The difference between petot and (\*p).peTot() is 1, which exceeds 1e-9, where
petot evaluates to 1,
(\*p).peTot() evaluates to 0, and
1e-9 evaluates to 1.0000000000000001e-09.
twoPatch metropolis
twoPatch wltmmc

There are some basic analysis tools available in FEASST
#########################################################
For reweighting, you can use:

[.../tools]$ ./rw.py --help
usage: rw.py [-h] [--inFile INFILE] [--outFile OUTFILE] [--lnz LNZ]
             [--phaseBoundary PHASEBOUNDARY]

optional arguments:
  -h, --help            show this help message and exit
  --inFile INFILE, -i INFILE
                        input collection matrix file
  --outFile OUTFILE, -o OUTFILE
                        output collection matrix file
  --lnz LNZ, -z LNZ     ln(activity)
  --phaseBoundary PHASEBOUNDARY, -p PHASEBOUNDARY
                        assign number of molecules as phase boundary

where the input is a collection matrix file. If an activity is not specified then it attempts to find two peaks for reweighting to phase equilibria. You can manually set the 'phase boundary' or let it automatically attempt to find the minimum between the two peaks. There can be issues if the lnPi is not well converged as has many local min/max, and its likely that some of Nate's python scripts on github are more sophisticated for these kind of special cases, minimum finding, etc. The tool assumes there is a file 'tmp/rstspace' to instantiate the Space object but I think all it wants to know is the number of particles or maybe the volume so it can output coexistence densities, so you could probably make an empty space object and still do the reweighting just fine (and its a simple calculation to check).

Debugging with GDB
####################

gdb is a very useful debugging tool, especially for identifying segfaults via backtraces. The -g flag in compliation pulls the symbols so that you can get correct line numbers in the gdb output.

In bash

.. code-block:: bash

   gdb [program executable name]
   r [flags]

gdb can also be used with python as

.. code-block:: bash

   gdb python
   r [python script] [optional flags]


TODO LIST
#####################

* json reader to server as 'input script' to launch simulations
* json as checkpoint file
* MD with stochastic dynamics integrator
* Perfect checkpointing
* Automated full-checkpoint testing
* remove printPressure from mc/criteria, printBeta, pairOrder, floppyBox, etc
* on the fly WL/TM lnPI error analysis ... accumulate 3 lnPIs by spliting each trial to each individual criteria class. Use them to compute all sorts of quantities.
* for xyz2bin, in afterAttempt MC, use unique hash on log file and xyz configuration for error check
  -- implmement with WLTMMC, use criteria to find order param column in log, then readxyz hash, find log hash match, demix conf based on the bin
* have criteria class backup colmat/stats periodically, based on sweeps?, that can be post processed (e.g., energy stats)
* combine pair_square_well, pair_hs, pair_hard_circle
* remove periodicity from x/y/z dimensions (no rush here)
* split functions.h into a variety of base_fileio, base_math, base_utils, etc
* pairhybrid rCut should be taken from pairVec, or atleast rCutMaxAll
* Eventually convert all raw pointers to shared pointers, which also allows removal of space from MC class
* PairHybrid also doens't need a space pointer.
* Use Histogram class for CriteriaWLTMMC instead of its own hard-coded version
* To reduce the size of Space, have it inherit multiple base classes, e.g.,
  Domain which contains box lengths and cell list, etc (but needs to know about particle positions?)
* Fix nomenclature.. atom == particle, mol == ?.. maybe change to sites / particles
* add ASSERT(rCutij.size() == 0 for linearShift in PairLJMulti so people don't run into issues with rCutij.clear
* Numerical implementation of quadratic equation coudl help with config bias: https://en.wikipedia.org/wiki/Quadratic_equation#Quadratic_formula_and_its_derivation
* Improve handling of default parameters for documentation and perhaps json (e.g. checkpointing above)?
* Move Add/mod new classes to API with links from README to API
* Combine PairLJCoulEwald and PairLJCoul in some way which doesn't involve copied code?
* change initEnergy in most implementations to use Inner() and reduce code complexity/copied code.
* nightly build -> unittests, test cases, coverage, valgrind, profiling, docs, python
* implement arbitrary order parameters as a class/factory method within CriteriaWLTMMC to allow users to define their own order parameters. These order parameters also must operate on Space/Pair objects (and also perhaps a Trial for expanded ensemble).
* runNumSweeps instead should have something where one generates the clones as vector<shrptr>, then runNumSweep takes these as input. That way one can modify the clones as one sees fit (also in multiprocessor restarts) before running the clones. It would take a lot of the hidden magic out without complicating the interface too drastically.
* move xdrfile and others to extern, change location of xdrfile file library away from "home" directory
* make extern/README.rst and others part of the documentation.
* Fix GSL memory leaks
* I prefer segfault on error for backtrace, but I should make it so all my packages do not segfault on back trace (CMakeList.txt macro?)

