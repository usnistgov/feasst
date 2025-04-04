==========
Tutorials
==========

The quickest way to get started with FEASST is to find a tutorial that is closest to what you would like to accomplish.
After verifying that you can reproduce the expected result, use the :doc:`../plugin/text_interface` documentation to better understand and modify the classes and arguments.
In most tutorials, a Python script is used to generate a text file based on input arguments to a formatted string.
The generated text file will be easier to understand than the Python script, and the text files can be run directly as `/path/to/feasst/build/bin/fst < input.txt`.
Include the text file when you :doc:`../CONTACT` us with issues.

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
   ../plugin/prefetch/tutorial/tutorial*
   ../plugin/aniso/tutorial/tutorial*
   ../plugin/fftw/tutorial/tutorial*
   ../plugin/morph/tutorial/tutorial*
   ../plugin/gibbs/tutorial/tutorial*
   ../plugin/netcdf/tutorial/tutorial*

Text file interface
======================

The above tutorials feature the :doc:`../plugin/text_interface` of FEASST, which is the recommended interface for most users.
The :doc:`../plugin/text_interface` comprises the names of many classes and their constructor arguments.
Thus, most users will want to search the :doc:`../plugin/text_interface` documentation for these classes.

Python and C++ interface
===============================

FEASST may also be called directly in Python and C++.
The :doc:`../python/README` and :doc:`../plugin/README` are not recommended for most users, and may change with minor version.

