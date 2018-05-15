SPC/E: muVT WL-TMMC
*******************************************************************************************************

.. code-block:: bash

    AUTO_GEN_DIR

Compute the equation of state and coexistence of SPC/E water in the grand canonical ensemble using flat-histogram methods.
Reproduce the results found in the NIST SRSW:

http://www.nist.gov/mml/csd/informatics_research/srsw.cfm

https://www.nist.gov/mml/csd/informatics/sat-tmmc-liquid-vapor-coexistence-properties-spce-water-lrc

The SPC/E model and reference configurations are provided at

https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10%C3%A5-cutoff

* `<testcase/spce/srsw/eostmmc/t525/run_*>`_ are example scripts to compile and run muvttmmcspce.
* ``h2o.525K*`` contain results from the SRSW.
* ``colMat`` file with the macrostate probability distribution (:math:`ln\Pi`), potential energy and collection matrix

To compare with the NIST SRSW:

* colMat columns 1:2 with ``h2o.525Kr*.p_macro.dat``

Now, in order to compare coexistence properties you will need to make use of histogram reweighting, as described in :doc:`1_lj_3_eostmmc_README`.
