*************************
Guidelines for Developers
*************************

Developers and contributors of FEASST may consider the following guidelines
in order to simplying collaboration.



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
* Consider feature compatability with swig python interface

Hotfix branch
**************

The hotfix branch is for implementation of bug fixes to the master branch.
These branches have a short life span. After mergeing with master and dev,
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

Find memory leaks with Valgrind
#################################

Valgrind helps to detect memory management bugs.

http://valgrind.org/

For example, to run valgrind on a particular test and output to text file

.. code-block:: bash

   valgrind ./unittest --gtest_filter=MC.* > out.txt 2>&1

* For unitilized value errors, try --track-origins=yes
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

Clean up with cpplint and pylint
#######################################

https://google.github.io/styleguide/cppguide.html

https://github.com/google/styleguide/tree/gh-pages/cpplint

Philosophical questions
#######################################

Why are there so many comments in the header files?
****************************************************

Because I assume most users are primarily concerned with the interface.
I prefer that one can understand most aspects of the class simply by reading the header file.
If you would like to contribute code, then document how you are most comfortable.

Why use camel case instead of underscores?
*******************************************

It saves one key stroke versus a "_" and one horizontal space.
Perhaps for the sake of readability, FEASST will switch to a more pythonic style but that would take a lot of work.
Feel free to stick to your favorite style in your contributed work, although ideally there is only one style per file!




