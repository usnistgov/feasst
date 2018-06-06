LJ: muVT TMMC
*******************************************************************************************************

.. code-block:: bash

    AUTO_GEN_DIR

Compute the equation of state and coexistence of Lennard-Jones particles in the grand canonical ensemble.
Reproduce the results found in the NIST SRSW:

http://www.nist.gov/mml/csd/informatics_research/srsw.cfm

http://mmlapps.nist.gov/srs/LJ_PURE/eostmmc.htm

* `<testcase/lj/srsw/eostmmc/1.2/run_*>`_ are example scripts to compile and run muvttmmclj.
* ``lj.msdb.*`` contain the results from the SRSW.
* ``colMat`` file with the macrostate probability distribution (:math:`ln\Pi`), potential energy and collection matrix

To compare the results with the NIST SRSW:

* colMat columns 1:2 with ``lj.msdb.t120.*.p_macro.dat``
* colMat columns 1:3 with ``lj.msdb.t120.*.energy.dat``

Now, in order to compare coexistence properties of the LJ potential against the NIST SRSW:

https://www.nist.gov/mml/csd/chemical-informatics-research-group/sat-tmmc-liquid-vapor-coexistence-properties-long-range

You will need to make use of histogram reweighting.

Histogram reweighting
=====================

For microcanonical partition function, :math:`\Omega(N,V,U)`, canonical, :math:`Q(N,V,T)` and grand canonical, :math:`\Xi(\mu, V, T)`, the probability to observe a state, :math:`\Pi` is given by

.. math::

   \Pi(U) = \frac{\Omega(N,V,U)e^{-\beta U}}{Q(N,V,T)}
 
   ln\Pi(U) + \beta U = ln\Omega(N,V,U) + Const
 
   \Pi(N, U) = \frac{\Omega(N,V,U)e^{-\beta U + \beta \mu N}}{\Xi(\mu,V,T)}
 
   ln\Pi(N,U) + \beta U - \beta \mu N = ln\Omega(N,U,V) + Const
 
   ln\Pi_{expected}(N, U; \mu', V, T') = ln\Pi(N, U; \mu, V, T) -\beta` U + \beta U + (\beta` \mu` - \beta \mu)N + Const
 
   lnz(\mu, \beta) = \beta\mu - ln(Lambda(\beta)^3)
 
   ln\Pi_{expected}(N; \mu', V, T) = ln\Pi(N; \mu, V, T) + N(lnz` - lnz) + Const

For histogram reweighting, apply the equation above to the macrostate probability distribution, :math:`ln\Pi`, and renormalize, :math:`\sum \Pi = 1`.
When there are multiple peaks in :math:`ln\Pi` that do not change greatly with system size, this often corresponds to different macroscopic phases.
One chooses the minimum in :math:`ln\Pi` as the boundary between these two phases.
The two phases are at coexistence when the sum of their probabilities are equal.
Extensive properties, :math:`A` are obtained first as a series of canonical ensemble averages for each value of :math:`N`, and then post-processed as a function of :math:`lnz` as grand canonical ensemble averages, :math:`A(\mu,V,T) = \sum_{N=0}^{N_{max}} A(N,V,T) \Pi(N; \mu,V,T)`

Example histogram reweighting tools are provided in `<tools/rw.cc>`_ and `<tools/rw.py>`_.
For C++, ``cd tools`` ``./compile.sh rw`` to make the ``rw`` executable.

To reweight :math:`ln\Pi` to different values of :math:`lnz`, use ``tools/rw -i colMat -z [lnz value]`` or ``tools/run.sh tools/rw.py -i colMat -z [lnz value]`` and the result is in ``colMatrw.txtrw``.
If your :math:`ln\Pi` has multiple peaks, you can obtain saturation by ``tools/rw -i colMat`` or ``tools/run.sh tools/rw.py -i colMat`` and the resulting coexistence properties are printed and the reweighted :math:`lnz` is in ``colMatrw``.

What chemical potential (lnz) should I use?
===========================================

It almost doesn't matter because you can reweight and the flat-histogram method gets over the barriers.
But you want to pick one that is closest to the properties of interest.
For example, if you want to study vapor-liquid coexistence then the closer you are to the actual coexistence lnz, the better.
Running short simulations will give you an idea if your simulation is favoring vapor too heavily (large peak for small N and monotonically decreasing :math:`ln\Pi`) or favoring liquid too heavily (large peak for large N and monotonically increasing :math:`ln\Pi`).
Plus, you can reweight your short simulation to see which lnz will get you closer to the region of interest.

What maximum number of particles should I use?
=================================================

The maximum number of particles is tricky.
You are essentially choosing the truncation of the grand canonical average sum provided in the `Histogram reweighting`_ section above.
As N gets large, :math:`ln\Pi` gets very small and doesn't significantly contribute to the grand canonical averages.
Thus, a simple rule of thumb is to say that once your :math:`ln\Pi` is about 20 units lower than the "peak" of interest in the macrostate distribution, then you no longer need values of N beyond this limit.
But you also have to remember than :math:`ln\Pi` changes when you are reweighting.
While it is more conservative to choose larger values of the maximum number of particles, it can also significantly increase the computational cost of the simulation.
Breaking the macrostate distribution into smaller pieces to perform the simulations helps, but you also want to make sure that your high density states are not simply stuck at some metastable state.

