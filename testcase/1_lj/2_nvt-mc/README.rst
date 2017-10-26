LJ: NVT MC
**************************************************************************************

.. code-block:: bash

    AUTO_GEN_DIR

Compute the average potential energy of Lennard-Jones particles in the canonical ensemble.
Reproduce the results found in the NIST SRSW:

https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm

In this case, implementation of a custom Analysis class to compute the average potential energy is shown.

Note that while this example is in the spirit of a monkey patch in python, implementing
a monkey patch on the SWIG python objects requires editing the vtable.
In this case, it may be easier to add the custom analysis in the source directory.

