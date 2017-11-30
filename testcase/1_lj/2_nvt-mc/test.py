"""
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http:#pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
"""

import unittest
import feasst

class TestLJ_SRSW_NVTMC(unittest.TestCase):
    def test(self):
        space = feasst.Space(3)
        rho = 1e-3  # number density
        nMol = 500     # number of particles
        space.initBoxLength((float(nMol)/rho)**(1./3.))   # set the cubic PBCs
        molNameSS = space.install_dir() + "/forcefield/data.lj"
        space.addMolInit(molNameSS)
        pair = feasst.PairLJ(space, 3,  # potential truncation at 3
            feasst.args({"cutType" : "lrc"}))
        pair.initEnergy()
        temperature = 0.9
        criteria = feasst.CriteriaMetropolis(1./temperature, 1.)
        mc = feasst.MC(space, pair, criteria)
        feasst.transformTrial(mc, "translate", 0.1)
        mc.nMolSeek(nMol)
        mc.initLog("log", int(1e4))
        mc.initMovie("movie", int(1e4))
        mc.initRestart("tmp/rst", int(1e4))
        mc.setNFreqTune(int(1e4))
        mc.runNumTrials(int(1e7))   # run equilibration

        # Run the production simulation and compute statistics on potential energy
        mc.setNFreqTune(0)  # do not tune during production
        pe = feasst.Accumulator()
        from itertools import islice, count
        for itrial in islice(count(1), int(1e7) - 1):
            mc.runNumTrials(1)
            pe.accumulate(pair.peTot()/float(space.nMol()))

        # Check average energy against the NIST SRSW
        # https:#mmlapps.nist.gov/srs/LJ_PURE/mc.htm
        # https:#www.nist.gov/programs-projects/nist-standard-reference-simulation-website
        peAv = pe.average()
        peStd = pe.blockStdev()
        peSRSW = -9.9165E-03
        peSRSWstd = 1.89E-05
        self.assertAlmostEqual(peAv, peSRSW, delta=2.576*(peSRSWstd+peStd))

if __name__ == "__main__":
    unittest.main()
