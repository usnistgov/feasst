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

class TestLJ_SRSW_REFCONF(unittest.TestCase):
  def test(self):
    space = feasst.Space(3)
    pair = feasst.PairFermiJagla(space, feasst.args({"rCut" : "3"}))
    molNameA = space.install_dir() + "/forcefield/data.lj"
    pair.initData(molNameA)
    molNameB = space.install_dir() + "/forcefield/data.ljb"
    pair.initData(molNameB)

    # set the cut-off
    pair.equateRcutForAllTypes()

    # set the pairwise B_0 parameter
    pair.epsijset(0, 0, 1.)
    pair.epsijset(0, 1, 1.)
    pair.epsijset(1, 1, 1.)

    # create clones of Space and Pair to perform two separate tests
    space2 = space.clone()
    pair2 = pair.clone(space2)

    # first, test the interaction between two particles
    xAdd = feasst.DoubleVector(space.dimen())
    pair.addMol(xAdd, molNameA)
    r = 1.0245
    xAdd[0] = r
    pair.addMol(xAdd, molNameB)
    pair.initEnergy()
    self.assertAlmostEqual(pair.peTot(), -0.24869232467128419, 15)

    # second, test a reference configuration
    pair2.epsijset(0, 0, 0.)
    pair2.epsijset(0, 1, 1.3218525572535642)
    pair2.epsijset(1, 1, 0.)
    config = space.install_dir() + \
      "/tutorial/2_jagla/1_ref-config_fermi/test_case.xyz"
    conf_file = feasst.make_ifstream(config)
    nA = 162
    for n in range(nA): pair2.addMol(molNameA)
    for n in range(nA): pair2.addMol(molNameB)
    pair2.readXYZ(conf_file)
    pair2.initEnergy()
    self.assertAlmostEqual(pair2.peTot()/space2.nMol(), -1.8207714765944178, 1)
    pair2.epsijset(0, 0, 1.)
    pair2.epsijset(1, 1, 1.)
    pair2.initEnergy()
    self.assertAlmostEqual(pair2.peTot()/space2.nMol(), -2.0590968411085431, 15)

if __name__ == "__main__":
    unittest.main()
