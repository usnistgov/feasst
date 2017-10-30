"""
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http:#pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
"""

import feasst

space = feasst.Space(3)
rho = 1e-3  # number density
nMol = 500     # number of particles
space.lset((float(nMol)/rho)**(1./3.))   # set the cubic PBCs
pair = feasst.PairHardSphere(space)
pair.initData(space.install_dir() + "/forcefield/data.atom")
pair.initEnergy()
criteria = feasst.CriteriaMetropolis(1., 1.)
mc = feasst.MC(space, pair, criteria)
feasst.transformTrial(mc, "translate", 0.1)
mc.nMolSeek(nMol)
mc.initLog("log", int(1e4))
mc.initMovie("movie", int(1e4))
mc.initRestart("tmp/rst", int(1e4))
mc.setNFreqTune(int(1e4))
mc.runNumTrials(int(1e7))   # run equilibration

# Run the production simulation
mc.runNumTrials(int(1e7))

