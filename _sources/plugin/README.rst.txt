***************************
All Classes/plugins
***************************

FEASST simulations are conducted by creating a series of classes in a specific order.
The first word in each line of a FEASST input text file is the name of a class, followed by pairs of class arguments.
The classes listed here represent all the capabilities of FEASST.

A FEASST plugin is a collection of related classes.
Unnecessary plugins, especially those with onerous dependencies, can be removed from the FEASST_PLUGINS variable in /path/to/feasst/CMakeLists.txt.
Similarly, plugins can also be added this way.
In both cases, FEASST must be reinstalled.

See the example plugin as a template for creating your own class or plugin.

.. toctree::

   threads/README
   utils/README
   math/README
   configuration/README
   system/README
   monte_carlo/README
   models/README
   steppers/README
   flat_histogram/README
   patch/README
   mayer/README
   xtc/README
   chain/README
   shape/README
   confinement/README
   charge/README
   opt_lj/README
   cluster/README
   egce/README
   morph/README
   beta_expanded/README
   prefetch/README
   aniso/README
   example/README
   fftw/README
   netcdf/README
