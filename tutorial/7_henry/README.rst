Example Calculation of the Henry constant and infinite-dilution isosteric heat of adsorption of TraPPE Ethane in LTA All-silica Zeolite
*******************************************************************************************************************************************

.. code-block:: bash

    AUTO_GEN_DIR

This calculation is accomplished via Widom Insertion to obtain N=1 state spatial
averages, using FEASST to generate the ghost insertions and perform energy
calculations.
It is based on the Henry coefficient calculation described in "Toward Rational
Design of Metal–Organic Frameworks for Sensing Applications: Efficient
Calculation of Adsorption Characteristics in Zero Loading Regime" by L.
Sarkisov, <em>J Phys Chem C</em>, 116(4):3025-3033, 2012
(https://doi.org/10.1021/jp210633w)

Requirements: This FEASST Script uses a specific branch and commit of the FEASST
source code, release tag v0.6.0
(https://github.com/usnistgov/feasst/releases/tag/v0.6.0),
as it utilizes tutorial-specific functions.
Checkout the specific release of FEASST via:

git checkout v0.6.0

Authors: Daniel W. Siderius, PhD; Harold W. Hatch, PhD
