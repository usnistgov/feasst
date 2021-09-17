***************************
Developer guidelines
***************************

Developers and contributors of FEASST may consider the following guidelines in order to simplify collaboration.

Branch policies
=======================

Branches in Git are extremely useful for working on different features within the same code base and then merging them into a final product or release.

https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging

Using and naming branches in a way that is consistent may simplify this process.
The policies below are based loosely on the following blog.

http://nvie.com/posts/a-successful-git-branching-model/

Main branch
--------------------------------------------------------------------------------

This is the release branch.
This branch may only merge with develop and hotfix branches.

* `GTEST: Unittest C++`_ and all tests pass
* `Valgrind: Find memory leaks`_
* `GCOV and LCOV: Test coverage`_ which must be >85%
* `cpplint and pylint: Clean up styling`_
* SWIG and Python unittested
* `Document`_ new features
* Tag :code:`git tag -a <version number>` and :code:`git push --tags`

Develop branch
--------------------------------------------------------------------------------

This is the development branch.

Tested and complete implementation of new features are added here.

* Code must compile and gtests pass successfully
* Consider feature compatibility with swig python interface
* When merging into this branch, use --squash to add the entire feature with a single commit.
  Alternatively, you could also rebase your private branch.

Hotfix branch
--------------------------------------------------------------------------------

The hotfix branch is for implementation of bug fixes to the main branch.
After merging with main and develop, they are deleted.

Feature branches
--------------------------------------------------------------------------------

These branches are for the development of new features.
These branch names must begin with characters which identify the main developer of the feature (e..g, hwh/feature)
The main developer sets the branch policy.

Dead (dead_*) branches
--------------------------------------------------------------------------------

These branches may record an incomplete attempt of a feature that may be relevant in the future.
These branch names must begin with the characters "dead".

* No rules. Code may not compile.
  HINT: rename branches with "git branch -m <newname>"

Pull requests
--------------------------------------------------------------------------------

To create a pull request, fork the usnistgov repo, create a new branch with your changes, and add the pull request, as detailed in https://opensource.com/article/19/7/create-pull-request-github

To incorporate the pull request into feasst
- git fetch usnistgov pull/ID/head:BRANCHNAME
- git checkout BRANCHNAME
- [make local changes]
- git push usnistgov BRANCHNAME

Tools
================================================================================

GTEST: Unittest C++
--------------------------------------------------------------------------------

The main branch requires GTEST coverage for all cpp plugins in plugin/name/test.

.. code-block:: bash

    cmake -DUSE_GTEST=ON ..
    make unittest -j12
    ./bin/unittest

* use :code:`--gtest_filter=*Name*` to run specific tests
* use :code:`./bin/unittest \|& more` to pipe stderr to stdout.
* use :code:`--gtest_shuffle` to randomize the order of the tests
* use :code:`--gtest_random_seed=SEED` to reproduce an specific order.

GDB or LLDB: Debugging
--------------------------------------------------------------------------------

gdb (or lldb on macOS) is especially useful for identifying segfaults via backtraces.
The -g flag in compilation pulls the symbols so that you can get correct line numbers in the gdb output.
Often, the optimization flags (e.g., -O3) can obfuscate the backtrace.
If that is the case, recompile without optimization (using cmake).

In bash

.. code-block:: bash

   gdb [program executable name]
   r [flags]

or

.. code-block:: bash

   gdb --batch --command=../dev/test.gdb ./bin/unittest


gdb can also be used with python as

.. code-block:: bash

   export PYTHONPATH=$PYTHONPATH:~/feasst/build/
   gdb python
   r [python script] [optional flags]

* use 'gdb catch throw' or 'lldb break set -E C++' to backtrace exceptions

* use gdb as a profiler by ctrl c in the middle and backtrace: https://stackoverflow.com/a/378024
* use gdb as a parallel profiler: http://poormansprofiler.org/

Valgrind: Find memory leaks
--------------------------------------------------------------------------------

Valgrind helps to detect memory management bugs.

http://valgrind.org/

For example, to run Valgrind on a particular test and output to text file

.. code-block:: bash

   valgrind ./unittest --gtest_filter=MC.* > out.txt 2>&1

* For uninitialized value errors, try --track-origins=yes
* For leaks, try --leak-check=full --show-leak-kinds=all
* Don't use profiler for leak checks. OMP causes "leaks" O.K.
* For suppress false-positives (e.g., gomp or gsl), use --gen-suppressions=all to generate suppression files

GCOV and LCOV: Test coverage
--------------------------------------------------------------------------------

GCC compilers allow testing of coverage with gcov and lcov for visualization.

* Code: currently implemented with Travis CI and CodeCov and available online.
  See .travis.yml for example of how to use lcov
* Use GCOV with CMake: cmake -DUSE_GCOV .
  Note: this disables optimization, so don't use it for production simulations.
* make coverage
* Open coverage/index.html in your browser.
* Go into "src" and ignore the external library coverage.

CCACHE: Speed up compilation time
--------------------------------------------------------------------------------

Something as trivial as changing a comment in a header file can lead to a massive recompile of the entire source.
Your previous compile is remembered by ccache, leading to near instant recompilation in the above example.

cpplint and pylint: Clean up styling
--------------------------------------------------------------------------------

https://google.github.io/styleguide/cppguide.html

https://github.com/google/styleguide/tree/gh-pages/cpplint

Document
================================================================================

Setup
--------------------------------------------------------------------------------

pip install sphinx breathe
doxygen with GENERATE_XML
run sphinx-quickstart, enable autodoc
add something like the following to your sphinx index.rst::

    .. doxygenclass:: Nutshell
       :project: nutshell
       :members:

add the following to your sphinx conf.py
  extensions = [ "breathe", "nbsphinx" ]
  breathe_projects = {"FEASST":"../xml"}
  breathe_domain_by_extension = {"h" : "cc"}

pip install sphinx_rtd_theme nbsphinx

run sphinx: make html

apt install graphviz graphviz-dev
pip install pygraphviz

Sphinx/Breathe/Doxygen notes
--------------------------------------------------------------------------------

* Link from rst file to C++ function: ``:cpp:func:`link <feasst::className::function()>```
* Link from rst file to C++ class: ``:cpp:class:`link <feasst::className>```
* Link from rst file to fst file: ``:doc:`/tutorial/asdf``` [note, / references root]
* Link from rst file to ipynb file : ```Tutorial <tutorial/tutorial.html>`_``
* Link from C++ to C++: ``className::function()``
* Link from C++ source to rst file: ``<a href="tutorial/asdf.html">test</a>``
* For math in C++ comments::

   \f$ latex code here \f$

* For tables, see monte_carlo/include/trial_compute_add.h

Pip notes
-------------------------

dev/tools/pip_install.sh

Style
================================================================================

Reference guides for C++
--------------------------------------------------------------------------------

* http://www.cplusplus.com/
* https://google.github.io/styleguide/cppguide.html
* http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines

Naming
--------------------------------------------------------------------------------

* ClassNames are mixed case with starting upper case letter
* member_names are lower case with underscores
* private_member_names\_ end with an underscore
* function_names are also lower case with underscores
* bools syntax: is_[accepted.., etc]
* MACROS and CONSTANTS are all upper case.
* Avoid MACROS and CONSTANTS.
* use "and", "or" instead of "&&", "||" (HWH: change this to follow Google?)

Functions
--------------------------------------------------------------------------------

* Use return values. Argument ordering: input (value or constant reference), then output (pointer only)
* Overloaded functions -> can you document all in a single comment? good
* No Default parameters on virtual functions

Classes
--------------------------------------------------------------------------------

* Nearly all data members should be private. Limit protected members
* member_name() returns const member
* set_member_name(member_name) sets member
* For setters with multiple arguments, the first are vector indices as in order x[0] = 3...
* getptr_member_name returns constant pointer (optimization only)

Loops and if
--------------------------------------------------------------------------------

* use of "for (auto element : container) { ... }" is dangerous
* for simple loops over containers, use "for (element : container)"
* for loops where you need the index, use:
  for (int index = 0; index < static_cast<int>(container.size()); ++index)

Auto
--------------------------------------------------------------------------------

* only use auto when the type is clear such as auto var = std::make_shared<..>.

Arguments
--------------------------------------------------------------------------------

* All arguments are provided as strings and converted to the expected type.
* Check that all arguments are used (e.g., like implicit none, a typo is caught).
* Argument defaults need to be set and clearly commented.
* If no default, it is a required argument.

Serialization
--------------------------------------------------------------------------------

* guided by https://isocpp.org/wiki/faq/serialization
* For inheritance hierarchy, a static deserialize_map is used to relate class
  name to template.
* Each object serializes a version that can be used for checks and backwards
  compatibility.
* utils_io.h contains many function templates for serialization.
* In particular, feasst_deserialize_fstdr() needs to be fixed.
* Don't forget to serialize (private) member data in new implementations.
* To compare differences between two serializations, paste into file and using "s/ /\r/g"

File output
--------------------------------------------------------------------------------

* comma-separated values (CSV) are the preferred format (e.g., comma deliminter)

For quick reference
================================================================================

* line counts [find . -name '*.cpp' -o -name '*.h' | xargs wc -l | sort -n]
* tutorial errors [ find . -name 'tutorial_failures.txt' | xargs cat ]
* tutorial errors [ for fl in `find . -name 'tutorial_failures.txt'`; do echo $fl; cat $fl; done ]
* clear tutorial errors [ for fl in `find . -name 'tutorial_failures.txt'`; do echo $fl; rm $fl; done ]
* clean docs before running depend.py again [ for dir in `ls --color=never -d *`; do rm $dir/doc/*rst; done ]
* find difference in serialization string: [ diff -u f1 f2 |colordiff  | perl /usr/share/doc/git/contrib/diff-highlight/diff-highlight | more ]

To Do List
================================================================================

* implement gibbs ensemble
* Make utils:lj,spce,etc derived classes of System ?
* Create benchmarking profile to compare among versions
* fix dependency linkers required by clang/cmake on macOS but not g++ on ubuntu
* consider optimization of Ewald: init ewald storage on particle types, precompute property index.
* when selecting from cpdf, use lnp instead of p?
* Python debug script: easy for user to run gdb/valgrind
* Toggle more debug levels, and localized to certain files/plugins, etc
* force precompute when reinitializing system, criteria, etc in MonteCarlo
* MonteCarlo subclass Simulation
* swig python wrap arguments std::istream and std::ostream for serialization
* add citations to tutorials (reweighting, etc) and also citation suggestions for MC objects
* VisitModels may prefer to update select properties (e.g., cell, eik)
* Jupyter notebook output should go to cells, not terminal that runs jupyter.
* lint file_[xyz,lmp]
* regrow but within near existing, for 'free dof, e.g. azimuthal in  angle, sphere in bond, etc'
* put cell list in finalize-heavy paradigm, update_positions updates cell of selection, finalize updates entire cell list. linked list
* config could use revert,finalize to update cell list only on finalization, and maybe not have to exclude from cell properties (why exclude?). same with ewald
* Rename TrialSelect->SelectTrial, TrialCompute->ComputeTrial. Rename Compute->Decide?.
* Somehow, trial_growth_expanded.h doesn't include debug.h but can compile with ASSERT
* Speed up RNG by maintaining int_distribution like dis_double
* Document utils lj, spce, rpm in tutorials
* Add a FAQ for common compile errors: "no known conversion from brace-enclosed initializer list to argtype" often means that a parameter was not converted to a string.
* Make a CachedRandom and CachedPotential for prefetch and avoid if statements that could slow down serial simulations.
* remove tutorial/fh.py
* Analyze/ModifyFactory optimization: use steps_per in factory to limit number of checks
* implement timer for profiles (with hierarchies by class... tried this, but its too slow. Time only infrequently?)
* implement a timer to auto-balance trial weights based on cpu time.
* More documentation/tutorial on how to create your own plugins, classes, etc
* add orientation argument to shapes with internal coordinate transformation
* System should track current energy of every potential for analysis (Criteria running energies may contain a part from each potential to simplify debugging).
* Consider using new state instead of old state in acceptance derivations
* Sort selection_of_all, or impose sorting in Select::add_particles. Currently, this leads to issues.
* Rename xyz files, and/or document more cleary (second line in xyz).
* Rename plugin chain->config_bias ?
* in optimizing where config only updates when trial finalized, how to build off new perturbed config in CB?
* EnergyMapNeighborCriteria should be able to take index of energy map stored by System
* Optimize TrialRemove for new_only by not computing interactions with neighbors
* Rename ewald plugin to charge
* Tunable implementation of configurational bias. When param is 0, rebuilds/renormalizes particles to prevent drift in bond lengths/angles.
* Add more documentation for Trials (TrialStage: Perturb/Select, TrialCompute. When finialize? reference?)
* Early tutorials are in ipynb, but transition out of that into .py, include .py on the website
* When swig4 finally works, can use multimap for dictionaries with more than one of same first key
* Rename Movie->XYZ
* Rename Stepper. Rename steps_per to trials_per, so as not to be confused with CB steps. Or rename num_steps to num_cb_steps. reference_index to cb_reference.
