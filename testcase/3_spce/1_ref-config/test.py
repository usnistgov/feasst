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

class TestSPCE_SRSW_REFCONF(unittest.TestCase):
    def test(self):

        for srswConfig in range(1, 3):
            space = feasst.Space(3)
            space.initBoxLength(20.)
            pair = feasst.PairLJCoulEwald(space, 10.)
            pair.initData(space.install_dir() + "/forcefield/data.spce")
            pair.initKSpace(5.6,  # alpha*L
                            27)  # k^2 < k2max cutoff

            # Set atom-based cut-off for comparison to SRSW reference config
            # Warning: atom-based cut-offs will fail because an optimized routine for
            # molecule-based cutoff is used by Monte Carlo trials.
            pair.initAtomCut(1)
            for iType in range(space.nParticleTypes()):
                for jType in range(space.nParticleTypes()):
                    pair.rCutijset(iType, jType, pair.rCut())

            # Read the configuration
            conf_file = feasst.make_ifstream(space.install_dir() + \
                "/testcase/3_spce/1_ref-config/spce_sample_config_periodic" + \
                str(srswConfig) + ".xyz")
            pair.readxyz(conf_file)
            pair.erfTable().tableOff()  # turn off the error function table
            pair.initEnergy()

            # Compare with the SRSW
            KJperMol2Kelvin = 1e3/feasst.idealGasConstant
            if srswConfig == 1:
                self.assertEqual(100, space.nMol())
                self.assertAlmostEqual(pair.peLJ()*KJperMol2Kelvin, 99538.736236886805, delta=feasst.DTOL)
                self.assertAlmostEqual(pair.peLRC()*KJperMol2Kelvin, -823.71499511652178, delta=feasst.DTOL)
                self.assertAlmostEqual(pair.peQReal()*KJperMol2Kelvin, -558888.91972785105, delta=1e-9)
                self.assertAlmostEqual(pair.peQFrr()*KJperMol2Kelvin, 6270.0937839397402, delta=1e-11)
                # Note, this is the negative of the sum of E_self and E_intra in SRSW
                self.assertAlmostEqual(pair.peQFrrSelf()*KJperMol2Kelvin, 34699.373693030466, delta=1e-9)
                self.assertAlmostEqual(pair.peTot()*KJperMol2Kelvin, -488603.17839517148, delta=1e-9)
            elif srswConfig == 2:
                self.assertEqual(200, space.nMol())
                self.assertAlmostEqual(pair.peLJ()*KJperMol2Kelvin, 193712.42256683594, delta=feasst.DTOL)
                self.assertAlmostEqual(pair.peLRC()*KJperMol2Kelvin, -3294.8599804660735, delta=feasst.DTOL)
                self.assertAlmostEqual(pair.peQReal()*KJperMol2Kelvin, -1192948.6940245957, delta=1e-8)
                self.assertAlmostEqual(pair.peQFrr()*KJperMol2Kelvin, 6034.9503269097158, delta=1e-9)
                # Note, this is the negative of the sum of E_self and E_intra in SRSW
                self.assertAlmostEqual(pair.peQFrrSelf()*KJperMol2Kelvin, 69398.747386013289, delta=1e-9)
                self.assertAlmostEqual(pair.peTot()*KJperMol2Kelvin, -1065894.9284973294, delta=1e-8)

if __name__ == "__main__":
    unittest.main()
