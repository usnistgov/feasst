SPC/E: NVT MC
**************************************************************************************

.. code-block:: bash

    testcase/3_spce/2_nvt-mc

Compute the average potential energy of SPC/E water in the canonical ensemble.

Reproduce the results found in "Pressure calculation in polar and charged systems using Ewald summation: Results
for the extended simple point charge model of water" by Gerhard Hummer, Niels Gro/nbech-Jensen, and Martin Neumann.

https://doi.org/10.1063/1.476834

Note that there are some systematic deviations from the published result because of the follow differences:

* The published result used Ewald boundary conditions :math:`\epsilon=65`
* The wave vectors are specified as :math:`k^2 <= 38` instead of :math:`k^2 < 38`
