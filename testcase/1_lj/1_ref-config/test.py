"""
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
"""

import unittest
import feasst

class TestLJ_SRSW_REFCONF(unittest.TestCase):
    def test(self):
        space = feasst.Space(3)
        space.initBoxLength(8)
        pair = feasst.PairLJ(space, feasst.args(
            {"rCut" : "3",          # potential truncation at 3
             "cutType" : "lrc"}))   # long range corrections

        # create clones of Space and Pair to perform two separate tests
        space2 = space.clone()
        pair2 = pair.clone(space2)

        # first, test the interaction between two particles
        xAdd = feasst.DoubleVector(space.dimen())  # position to add is origin
        pair.addMol(xAdd)        # add first molecule
        r = 1.2345
        xAdd[0] = r              # position of second particle
        pair.addMol(xAdd)        # add second particle
        pair.initEnergy()        # compute energy
        peExact = 4*(pow(r, -12) - pow(r, -6))
        self.assertAlmostEqual(pair.peLJ(), peExact, 15)
        peExactLRC = (8./3.)*feasst.PI*space.nMol()**2/space.vol() \
          *((1./3.)*pair.rCut()**(-9) - pair.rCut()**(-3))
        self.assertAlmostEqual(pair.peLRC(), peExactLRC, 15)

        # second, compare with the reference configuration
        config = space.install_dir() + \
          "/testcase/1_lj/1_ref-config/lj_sample_config_periodic4.xyz"
        conf_file = feasst.make_ifstream(config)
        pair2.readXYZ(conf_file)                  # read the xyz file
        pair2.initEnergy()                        # compute energy
        peLJ = -16.790321304625856
        peLRC = -0.5451660014945704
        self.assertAlmostEqual(pair2.peLRC(), peLRC, 15)
        self.assertAlmostEqual(pair2.peLJ(), peLJ, 15)

if __name__ == "__main__":
    unittest.main()
