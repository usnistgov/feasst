***************************
Plugins
***************************

A plugin is a distinct piece of the FEASST simulation program.
They are designed to reduce dependencies and allow for others to easily develop their own plugins and modifications without affecting the rest of the code.
See the example plugin as a template for creating your own.
Plugins may be added or removed by changing the FEASST_PLUGINS variable in CMakeLists.txt of the root directory of FEASST.

.. toctree::

   utils/README
   math/README
   configuration/README
   system/README
   monte_carlo/README
   example/README
   models/README
   steppers/README
   flat_histogram/README
   patch/README
   mayer/README
   xtc/README
   chain/README
   shape/README
   confinement/README
   ewald/README
   opt_lj/README
   threads/README

Plugins in development
************************

.. toctree::

   cluster/README
   egce/README
   morph/README
   beta_expanded/README
   prefetch/README
