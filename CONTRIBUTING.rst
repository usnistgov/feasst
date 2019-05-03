***************************
Contributing
***************************

Developers and contributors of FEASST may consider the following guidelines in order to simplify collaboration.

Branch policies
=======================

Branches in Git are extremely useful for working on different features within the same code base and then merging them into a final product or release.

https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging

Using and naming branches in a way that is consistent may simplify this process.
The policies below are based loosely on the following blog.

http://nvie.com/posts/a-successful-git-branching-model/

Master branch
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

The hotfix branch is for implementation of bug fixes to the master branch.
After merging with master and develop, they are deleted.

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

Tools
================================================================================

GTEST: Unittest C++
--------------------------------------------------------------------------------

The master branch requires GTEST coverage for all cpp plugins in plugin/name/test.

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

In bash

.. code-block:: bash

   gdb [program executable name]
   r [flags]

gdb can also be used with python as

.. code-block:: bash

   export PYTHONPATH=$PYTHONPATH:~/feasst/build/
   gdb python
   r [python script] [optional flags]

* use 'gdb catch throw' or 'lldb break set -E C++' to backtrace exceptions

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

Sphinx/Breathe/Doxygen notes
--------------------------------------------------------------------------------

* Link from rst file to C++ function: ``:cpp:func:`className::function()```
* Link from rst file to fst file: ``:doc:`/tutorial/asdf```
* Link from C++ to C++: ``className::function()``
* Link from C++ source to rst file: ``<a href="tutorial/asdf.html">test</a>``
* For math in C++ comments::

   \f$ latex code here \f$

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

Auto
--------------------------------------------------------------------------------

* only use auto when the type is clear such as auto var = std::make_shared<..>.
* use of "for (auto element : container) { ... }" is dangerous

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

For quick reference
================================================================================

* line counts [find . -name '*.cpp' -o -name '*.h' | xargs wc -l | sort -n]

To Do List
================================================================================

* ideal gas as the first tutorial/testcase
* specify units in LMP data files?
* fix dependency linkers required by clang/cmake on macOS but not g++ on ubuntu
* detail branch rules (branch namespace and descriptive file, merge squash, release tags, etc).
* py/feasst.i depends on which plugins that you use. how to make this user friendly?
* consider a different way to interface selection and configuration.
* precompute long range corrections and faster types calculation
* multistage insertions and deletions
* neighbor lists
* windowing (MonteCarloFactory class)
* trial regrow to include grand canonical insertions

