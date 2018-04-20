*************************
Guidelines for developers
*************************

Developers and contributors of FEASST may consider the following guidelines
in order to simplify collaboration.

Semantic versioning
####################

FEASST version policy follows Semantic Versioning 2.0.0 http://semver.org/spec/v2.0.0.html.

In a nutshell, this means that version x.y.z-a-[version control id][abbreviated commit hash]

* x. Major version number. Public API is backwards incompatible with previous versions.
* y. Minor version number. New features but the public API is fully backwards compatible.
* z. This version number changes if the source code has changed in the slightest.
* a. This version changes if the documentation changes but the source code is the same.
* And the version control id is "g" for GitHub, for example, followed by the commit hash.

Although all efforts are made for restart files to be backwards compatible, they may only be compatible for the same version x.y.z.

Branch policies
###############

Branches in Git are extremely useful for working on different features within
the same code base and then merging them into a final product or release.

https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging

Using and naming branches in a way that is consistent may simplify this
process. The policies below are based loosely on the following blog.

http://nvie.com/posts/a-successful-git-branching-model/

Master branch
*************

This is the release branch.
This branch may only merge with development and hotfix branches.

* `Unittest with GTEST`_ and all tests pass
* `Find memory leaks with Valgrind`_
* `Test coverage with GCOV and LCOV`_ which must be >75% (aim for 100%)
* `Clean up with cpplint and pylint`_
* SWIG and Python unittested
* Update packages
* Document new features

Development (dev) branch
*************************

This is the development branch.

Tested and complete implementation of new features are added here.

* Code must compile and gtests pass successfully
* Consider feature compatibility with swig python interface

Hotfix branch
**************

The hotfix branch is for implementation of bug fixes to the master branch.
These branches have a short life span. After merging with master and dev,
they are deleted.

Feature (ftr_*) branches
*************************

These branches are for the development of new features.
These branch names must begin with the characters "ftr_*".
Once tested and complete, these branches may be merged into dev.

* Code must compile.
  HINT: use "git stash" over lazy commits, or use your own local branch

Dead (dead_*) branches
***********************

These branches may record an incomplete attempt of a feature that may be
relevant in the future.
These branch names must begin with the characters "dead".

* No rules. Code may not compile.
  HINT: rename branches with "git branch -m <newname>"

Public/private repositories with package tools
###############################################

.. include:: ../tools/package/README.rst

Unittest with GTEST
####################

The master branch requires GTEST coverage. All interface files have an
associated ``*_unittest`` file such as the following.

.. code-block:: c++

    #include "space.h"

    TEST(Space, init_config){
      Space space(3);
      EXPECT_EQ(3, space.dimen());  // 3D space
      EXPECT_FALSE(space.dimen() == 2);
      EXPECT_LT(2, space.dimen());
      EXPECT_NEAR(3., space.dimen(), DTOL);
    }

CMake and the Makefile will automatically detect and compile these tests.
For example

.. code-block:: bash

   cmake -DUSE_GTEST=ON .
   make unittest -j 8
   ./bin/unittest --gtest_shuffle --gtest_filter=Space.ini*

The optional gtest flags randomly shuffle the order of the tests, and run
specific tests, respectively.
To use a seed to reproduce an error, use --gtest_random_seed=SEED.

If there are a lot of errors, typically the best testing order would be::

    shape
    base
    functions
    table
    accumulator
    histogram
    random
    space
    pair
    analyze
    criteria
    trial
    mc
    ui_*

Debugging with GDB or LLDB
###########################

gdb (or lldb on macOS) is especially useful for identifying segfaults via backtraces. The -g flag in compilation pulls the symbols so that you can get correct line numbers in the gdb output.

In bash

.. code-block:: bash

   gdb [program executable name]
   r [flags]

gdb can also be used with python as

.. code-block:: bash

   gdb python
   r [python script] [optional flags]


Find memory leaks with Valgrind
#################################

Valgrind helps to detect memory management bugs.

http://valgrind.org/

For example, to run Valgrind on a particular test and output to text file

.. code-block:: bash

   valgrind ./unittest --gtest_filter=MC.* > out.txt 2>&1

* For uninitialized value errors, try --track-origins=yes
* For leaks, try --leak-check=full --show-leak-kinds=all
* Don't use profiler for leak checks. OMP causes "leaks" O.K.
* For suppress false-positives (e.g., gomp or gsl), use --gen-suppressions=all to generate suppression files

Test coverage with GCOV and LCOV
###################################

GCC compilers allow testing of coverage with gcov and lcov for visualization.

* Use GCOV with CMake: cmake -DUSE_GCOV .
  Note: this disables optimization, so don't use it for production simulations.
* make coverage
* Open coverage/index.html in your browser.
* Go into "src" and ignore the external library coverage.

Speed up compilation time with ccache
######################################

If you haven't used ccache before, give it a try.

Clean up with cpplint and pylint
#######################################

https://google.github.io/styleguide/cppguide.html

https://github.com/google/styleguide/tree/gh-pages/cpplint

Tag and archive
#####################

git tag -a <version number>

git archive master --prefix='feasst'`git describe master`'/' | gzip > feasst`git describe master`.tar.gz

Documentation setup
####################

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

Sphinx/Breathe/Doxygen tips
############################

* Link from rst file to C++ function: ``:cpp:func:`className::function()```
* Link from rst file to fst file: ``:doc:`/testcase/asdf```
* Link from C++ to C++: ``className::function()``
* Link from C++ source to rst file: ``<a href="testcase/asdf.html">test</a>``
* For math in C++ comments::

   \f$ latex code here \f$

To do list
#####################

* json reader to server as 'input script' to launch simulations
* json as checkpoint file
* MD with stochastic dynamics integrator
* Perfect checkpointing
* Automated full-checkpoint testing
* for xyz2bin, in afterAttempt MC, use unique hash on log file and xyz configuration for error check
  -- implement with WLTMMC, use criteria to find order param column in log, then readxyz hash, find log hash match, demix conf based on the bin
* have criteria class backup colmat/stats periodically, based on sweeps?, that can be post processed (e.g., energy stats)
* remove periodicity from x/y/z dimensions (no rush here)
* split functions.h into a variety of base_fileio, base_math, base_utils, etc
* Use Histogram class for CriteriaWLTMMC instead of its own hard-coded version
* To reduce the size of Space, have it inherit multiple base classes, e.g.,
  Domain which contains box lengths and cell list, etc (but needs to know about particle positions?)
* Fix nomenclature.. atom == particle, mol == ?.. maybe change to sites / particles
* add ASSERT(rCutij.size() == 0 for linearShift in PairLJMulti so people don't run into issues with rCutij.clear
* Numerical implementation of quadratic equation coudl help with config bias: https://en.wikipedia.org/wiki/Quadratic_equation#Quadratic_formula_and_its_derivation
* Improve handling of default parameters for documentation and perhaps json (e.g. checkpointing above)?
* implement arbitrary order parameters as a class/factory method within CriteriaWLTMMC to allow users to define their own order parameters. These order parameters also must operate on Space/Pair objects (and also perhaps a Trial for expanded ensemble).
* runNumSweeps instead should have something where one generates the clones as vector<shrptr>, then runNumSweep takes these as input. That way one can modify the clones as one sees fit (also in multiprocessor restarts) before running the clones. It would take a lot of the hidden magic out without complicating the interface too drastically.
* Fix GSL memory leaks
* add PairPatchKF loops to Pair base class for molecule-center-based cut-offs before loops through particles (as done for PairLJCoulEwald),
* Add verbosity level for printing debug messages
* Make a citation generator for techniques that are used.
* move checkBond from MC to analyze class for rigid particle simulations (or as part of checkE?)
* Move some analysis outside of space class
  - but what if things like trials use the analysis, such as clusters?
  - ""This trial requires MC to have an analyze class...""
* document log file and tuning, changing acceptance percentage
* Include PDF version of manual in release tarballs
* Package FEASST as a Docker app
* move unittests to testcases: ljmultilambda2D, NPTandLJvLJMulti, semigrand,
* Fix the API TOC for non-derived classes (e.g., TOC disappears when you select "Histogram")
* Add test case with AnalyzeCluster and AnalyzeScatter.
* Add test case with PairTabular1D (and/or improve beyond linear interpolation).
* Incorporate more "tools" as part of the test cases (e.g., rw, restart, xyz2bin, etc).
* feasst.nist.gov instead of pages.nist.gov/feasst ?
* sphinx_rtd_theme doesn't work with NIST header/footer
* Use short version number for display on FEASST html
* Customize restart directory name
* Fix testcase 1/3 by taking arguments with unittest in python.
* Use packages to dynamically add the disclaim in the release.
* rename feasst.str.. python doesn't like it.
* Eventually convert all raw pointers to shared pointers (but interface blind to shared_ptr vs unique_ptr), which also allows removal of space from MC class
* Custom derived Analysis class should be able to clone without more boiler plate
* Operations to assign atoms to groups, and be able to apply groups to most operations like analyze, similar to lammps
* Separate PairLJCoul into different classes for LJ, real-space Q and Ewald. Use PairHybrid to combine them.
* Once subspace efficiently implemented, pairHybrid can store different pairs with different subspaces to optimize e.g., the LJ interaction on oxygens only (for MFB?)
* Rename some variables... mpart should be msite? ran vs rand.
* Implement Mayer sampling within the context of CriteriaMayer + MC handling trials (e.g., need a different way to keep track of running energy due to overflow)


