***************************
Developer guidelines
***************************

Developers and contributors of FEASST may consider the following guidelines in order to simplify collaboration.

Branch policies
=======================

Branches in Git allow simultaneous work on different features within the same code base, and merging into a final product or release.
Consistent branch naming simplifies this process.

Main branch
--------------------------------------------------------------------------------

This is the release branch.
This branch merges with develop and hotfix branches.

* `GTEST: Unittest C++`_ and all tests pass
* `Valgrind: Find memory leaks`_
* `GCOV and LCOV: Test coverage`_
* `cpplint and pylint: Clean up styling`_
* SWIG and Python unittested
* `Document`_ new features
* Tag :code:`git tag -a <version number>` and :code:`git push --tags`

Develop branch
--------------------------------------------------------------------------------

This is the development branch.
Tested and complete implementation of new features are added here.

* Code must compile and gtests pass successfully
* When merging into this branch, --squash to add the entire feature with a single commit.
  The feature branch can be saved by the developer as complete_feature_commit for future reference.
  Alternatively, you could rebase your private branch.

Hotfix branch
--------------------------------------------------------------------------------

The hotfix branch is for implementation of bug fixes to the main branch.

Feature branches
--------------------------------------------------------------------------------

These branches are for the development of new features.
These branch names begin with the main developer identity (e..g, username/feature)
The main developer sets the branch policy.

Dead (dead_*) branches
--------------------------------------------------------------------------------

These branches may record an incomplete attempt of a feature that may be relevant in the future.
These branch names must begin with the characters "dead".

* Code may not compile.
  Rename branches with "git branch -m <newname>"

Pull requests
--------------------------------------------------------------------------------

To create a pull request:

* fork the usnistgov repo
* create a new feature branch "user/feature" to implement your changes (git checkout develop; git checkout -b user/feature)
* make as many commits as you want for yourself, which will be squashed into one commit when you're finished as follows
* when your feature is ready and well tested and you're ready to submit a pull request, first you can check if any more changes were made to develop (git checkout develop; git pull usnistgov develop; git checkout user/feature; git merge develop)
* squash your commits into a new branch "testmerge" off develop (git checkout develop; git checkout -b testmerge; git merge --squash user/feature)
* submit the "testmerge" branch in the pull request.
* You can keep your "user/feature" branch with more fine-grain commits for your own purposes or later testing if an issue is discovered. I often rename these "complete_feature_[commit]" after merging to develop where the commit is 6-10 of the first characters of the squash merge commit (for later record).

Try to copy the style of existing commits in the main branch. This means reducing the number of commits in a pull request (e.g., consider a git merge), and describe commits in a similar style as existing commits. For example, the description should be a single sentence that is general and not redundant. The commit description should be concise (e.g., the classes or plugins changed or added).

For lead developers to incorporate the pull request into feasst
- git fetch usnistgov pull/ID/head:BRANCHNAME
- git checkout BRANCHNAME
- [make local changes if required. Can "git commit --amend"]
- If authorship is missing after a squash merge, use git commit --amend --author="My Nick my.adress@email.com"
- git push usnistgov BRANCHNAME

Tools
================================================================================

GTEST: Unittest C++
--------------------------------------------------------------------------------

The main branch requires GTEST coverage for all cpp plugins in "plugin/[name]/test."

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
The -g flag in compilation is required for correct line numbers in the gdb output.
Often, the optimization flags (e.g., -O3) can obfuscate the backtrace.
If that is the case, recompile without optimization (using CMakeLists.txt).

In BASH

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

* backtrace exceptions with 'gdb catch throw' or 'lldb break set -E C++'
* profile by running a long simulation, pause in the middle with ctrl c, then backtrace

Valgrind: Find memory leaks
--------------------------------------------------------------------------------

Valgrind detects memory issues.
For example, to run Valgrind on a particular test and output to text file in BASH

.. code-block:: bash

   valgrind ./unittest --gtest_filter=MC.* > out.txt 2>&1

* For uninitialized value errors, try --track-origins=yes
* For leaks, try --leak-check=full --show-leak-kinds=all
* Don't use profiler for leak checks. OMP causes "leaks" O.K.
* Suppress false-positives (e.g., gomp or gsl) with --gen-suppressions=all to generate suppression files

GCOV and LCOV: Test coverage
--------------------------------------------------------------------------------

GCC compilers allow testing of coverage with gcov and lcov for visualization.

* Code: Travis CI and CodeCov and available online.
* Use GCOV with CMake: cmake -DUSE_GCOV -DUSE_GTEST=ON ..
  Note: this disables optimization, so don't use it for production simulations.
* make coverage
* Open coverage/index.html in your browser.
* Go into "src" and ignore the external library coverage.

CCACHE: Speed up compilation time
--------------------------------------------------------------------------------

Changing a comment in a header file can lead to a massive recompile of the entire source.
If your previous compile is ccache'd, recompilation after adding a comment is near instantaneous.

cpplint and pylint: Clean up styling
--------------------------------------------------------------------------------

https://google.github.io/styleguide/cppguide.html

https://github.com/cpplint/cpplint

Document
================================================================================

Setup
--------------------------------------------------------------------------------

sudo apt install doxygen pandoc
pip install sphinx breathe pandoc
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

pip install sphinx_rtd_theme nbsphinx sphinxcontrib-bibtex

run sphinx: make html

apt install graphviz graphviz-dev pandoc

pip install pygraphviz breathe pandoc

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

For Ubuntu 22, I had to comment out lines 713-714 of ~/.pyenv/feasst/lib/python3.10/site-packages/breathe/renderer/sphinxrenderer.py

                #assert isinstance(n, addnodes.desc_annotation)
                #assert n.astext()[-1] == " "

Sphinx bibtex
----------------

Add references in header file documentation as \ :footcite:p:`bibtex_name`\endrst or :footcite:t: to include author name.
and update the bibtex in /feasst/dev/sphinx/refs.bib

References in comment blocks must begin in the \rst environment and end with

References:

.. footbibliography::
\endrst

Pip notes
-------------------------

dev/tools/pip_install.sh

Style
================================================================================

Reference guides for C++
--------------------------------------------------------------------------------

* https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines
* https://google.github.io/styleguide/cppguide.html

Naming
--------------------------------------------------------------------------------

* ClassNames are mixed case with starting upper case letter
* member_names are lower case with underscores
* private_member_names\_ end with an underscore
* function_names are also lower case with underscores
* bools syntax: is_[accepted.., etc]
* MACROS and CONSTANTS are all upper case. (avoid MACROS and CONSTANTS).

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

* for simple loops over containers, use "for (element : container)"
* for loops where you need the index, use:
  for (int index = 0; index < static_cast<int>(container.size()); ++index)

Auto
--------------------------------------------------------------------------------

* use auto when the type is clear such as auto var = std::make_shared<..>.

Arguments
--------------------------------------------------------------------------------

* All arguments are provided as strings and converted to the expected type.
* Check that all arguments are used (e.g., like implicit none, a typo is caught).
* Argument defaults need to be set and clearly commented.
* If no default, it is a required argument.

Serialization
--------------------------------------------------------------------------------

* See https://isocpp.org/wiki/faq/serialization
* For inheritance hierarchy, a static deserialize_map is used to relate class name to template.
* Each object serializes a version that can be used for checks and backwards compatibility.
* utils_io.h contains many function templates for serialization.
* Serialize (private) member data
* To compare differences between two serializations, paste into a file and vim/sed "s/ /\r/g"

File output
--------------------------------------------------------------------------------

* comma-separated values (CSV) are preferred over space-separated.

For quick reference
================================================================================

* line counts [find . -name '*.cpp' -o -name '*.h' | xargs wc -l | sort -n]
* tutorial errors [ find . -name 'tutorial_failures.txt' | xargs cat ]
* tutorial errors [ for fl in `find . -name 'tutorial_failures.txt'`; do echo $fl; cat $fl; done ]
* launch errors [ for fl in `find . -name 'launch_failures.txt'`; do echo $fl; cat $fl | grep -v "Terminating because Checkpoint"; done ]
* clear tutorial errors [ for fl in `find . -name 'tutorial_failures.txt'`; do echo $fl; rm $fl; done ]
* clean docs before running depend.py again [ for dir in `ls --color=never -d *`; do rm $dir/doc/*rst; done ]
* screen html errors [ make html > tt 2&>1; grep -v "WARNING: document i" tt | grep -v "WARNING: Duplicate" | grep -v "Declaration is" > ttt ]
* find all headers in the public interface [ find . -name '*.h' | xargs grep "^  \/\*\* \@name Arguments$" ]
* find difference in serialization string: [ diff -u f1 f2 |colordiff  | perl /usr/share/doc/git/contrib/diff-highlight/diff-highlight | more ]

To Do List
================================================================================

* profile: Create benchmarking profile to compare among versions
* profile: implement timer for profiles (with hierarchies by class... tried this, but its too slow. Time only infrequently?)
* profile: implement a timer to auto-balance trial weights based on cpu time.
* debug: Toggle more debug levels, and localized to certain files/plugins, etc
* debug: Implement a class-specific debug output setting
* compile: fix dependency linkers required by clang/cmake on macOS but not g++ on ubuntu
* opt: when selecting from cpdf, use lnp instead of p?
* force precompute when reinitializing system, criteria, etc in MonteCarlo
* MonteCarlo subclass Simulation
* add citations to tutorials (reweighting, etc) and also citation suggestions for MC objects
* VisitModel may prefer to update select properties (e.g., cell, eik)
* Rename TrialSelect->SelectTrial, TrialCompute->ComputeTrial. Rename Compute->Decide?.
* Speed up RNG by maintaining int_distribution like dis_double
* add orientation argument to shapes with internal coordinate transformation
* Sort selection_of_all, or impose sorting in Select::add_particles. Currently, this leads to issues.
* Rename xyz files, and/or document more cleary (second line in xyz).
* Rename plugin chain->config_bias ?
* in optimizing where config only updates when trial finalized, how to build off new perturbed config in CB?
* Optimize TrialRemove for new_only by not computing interactions with neighbors
* Tunable implementation of configurational bias. When param is 0, rebuilds/renormalizes particles to prevent drift in bond lengths/angles.
* (repeat) regrow but within near existing, for 'free dof, e.g. azimuthal in  angle, sphere in bond, etc'
* Rename Movie->PrintXYZ
* Patch custom model params not present in mc.configuration().model_params (affects FileXYZPatch).
* early rejection scheme: https://doi.org/10.1080/00268976.2014.897392
* Windows with non-integer macrostates?
* For unknown reasons, VisitModelOuterCutoff had energy issues with RPM
* Add TrialParticlePivot to TrialGrow (randomly orients particle about site). Or, more generally, say num_steps=-1 combines stages into one.
* Support compressed trajectory formats: xtc, dcd, etc
* Implement Jeff's parallel method via CollectionMatrixSplice that allows exchange of window ranges with overlapping simulations
* Update Translate tunable maximum when volume changes
* Move Trial checks so that they can be applied to GhostTrialGrow
* Restart Prefetch. Does Run::num_trials work properly?
* optimize Ewald::update_wave_vector for NPT
* Reduce size of Checkpoint files for cell/neighbor lists (re-compute instead of checkpointing them). Also large tables.
* Represent relative rigid bodies as screw motion: https://en.wikipedia.org/wiki/Screw_theory
* Ewald mod k2max like LAMMPS
* Add precompute for BondFourBody, ThreeBody, etc, to speed up? But different dihedrals have different coefficients.
* Allow one script to contain multiple MonteCarlo, CollectionMatrixSplice, Prefetch in any order? (checkpointing is difficult)
* When trials start, check to see if there is a trial that uses weight_per_number_fraction but there are fixed particles (or, see if there are weight_per_number for all types unless excluded?)
* Add a FAQ for sim questions, such as, an overview of various table potential options, etc.
* Optimize BondVisitor that uses deserialize_map and strings in inner loop
* Add more documentation/examples of analyzing stdev of the mean with block analysis. Output individual block averages for custom analysis? Correlation time? Move Accumulator example to text interface. Expose Accumulator options (stepper takes Accumulator arguments).
* Have the tests override hours checkpoint, etc so that users don't have bad values
* Reduce input text file with class that factories the creation of analyze and modify (e.g., all using the same trials_per_iteration, file name prefixes, etc)
* fstprt files use label strings instead of numbers (document this, 0-O, 1-H for spce, etc)
* pip install feasst
* In class documentation, link to tutorials that use the class
* Every time trial is added, determine which molecules are excluded from weight_per_num_fraction
* make a gui/software/script that walks through the building of a FEASST input file.
* For 0.26, Remove Random::[time, default] arguments
* For 0.26, Remove ConvertToRefPotential, ProfileTrials, RemoveModify, RemoveAnalyze, RemoveTrial, Run::until_criteria_complete, Criteria/Stepper::iteration
* For 0.26, Depreciate and update AngleSquareWell::min/max to min_degrees/max_degrees
* For 0.26, Remove serialization version checks grep "if (version >"
