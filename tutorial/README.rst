==========
Tutorials
==========

The following tutorials are also located in the code repository with the same path as the html (e.g., in /path/to/feasst/plugin/[name]/tutorial/).

Basic tutorials:

.. toctree::
   :glob:

   tutorial.ipynb
   launch

Build and test models:

.. toctree::
   :glob:

   ../plugin/monte_carlo/tutorial/tutorial*
   ../plugin/charge/tutorial/tutorial*
   ../plugin/models/tutorial/tutorial*
   ../plugin/example/tutorial/tutorial*

Flat-histogram simulations:

.. toctree::
   :glob:

   ../plugin/flat_histogram/tutorial/tutorial*

Others:

.. toctree::
   :glob:

   ../plugin/confinement/tutorial/tutorial*
   ../plugin/chain/tutorial/tutorial*
   ../plugin/cluster/tutorial/tutorial*
   ../plugin/steppers/tutorial/tutorial*
   ../plugin/beta_expanded/tutorial/tutorial*
   ../plugin/egce/tutorial/tutorial*
   ../plugin/mayer/tutorial/tutorial*
   ../plugin/math/tutorial/tutorial*
   ../plugin/prefetch/tutorial/tutorial*
   ../plugin/aniso/tutorial/tutorial*
   ../plugin/fftw/tutorial/tutorial*
   ../plugin/morph/tutorial/tutorial*
   ../plugin/netcdf/tutorial/tutorial*

Text file interface
======================

The above tutorials feature the text-based interface of FEASST.
The text interface comprises the names of all Classes and their constructor arguments.
Thus, most users will want to search the documentation for these classes.
Documentation for a particular class can be found by searching the `list of all capabilities <../plugin/README.html>`_ (typically, using ctrl-f).
All available arguments for that class can then be found in the list of "args" of that class, or its base classes.

For developers, the bash command "grep" is also an option to search for more information on classes and their arguments.
For example, if you would like more information on `RandomMT19937 <plugin/math/doc/RandomMT19937.html>`_ but are not sure where to find it, you could search headers files

.. code-block:: bash

   grep -r --include=*.h RandomMT19937

And find that the class is part of the `Math <plugin/math/README.html>`_ plugin.

Python and C++ interface
===============================

FEASST may also be called directly in Python and C++ as a `library <library/tutorial.html>`_.
Note that the Python and C++ interfaces may change with minor version.
Thus, the Python and C++ interfaces are not recommended for beginning users.

