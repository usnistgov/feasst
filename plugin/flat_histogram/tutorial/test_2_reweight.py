import unittest
import pandas as pd
import reweight
import feasst

class TestFlatHistogramReweight(unittest.TestCase):
    """Test reweighting"""
    def test(self):
        """Compare with NIST SRSW
        https://www.nist.gov/programs-projects/nist-standard-reference-simulation-website
        https://www.nist.gov/mml/csd/chemical-informatics-research-group/sat-tmmc-liquid-vapor-coexistence-properties-long-range
        """
        dfsrsw = pd.read_csv("../test/data/stat120.csv")
        lnpi = dfsrsw["lnPI"]
        #print(dfsrsw)
        rw = reweight.reweight()
        rw.set_beta(1./1.2)

        # the statistical average slightly affected normalization
        self.assertAlmostEqual(1, rw.area(lnpi), 4)
        lnpi = rw.norm(lnpi)
        self.assertAlmostEqual(1, rw.area(lnpi), 15)

        # no phase boundary detected at this mu
        self.assertTrue(rw._nboundary(lnpi) is None)

        rw.set_lnpi(lnpi)
        mu = -2.902929
        mus = rw.lnzsat(lnpi, mu=mu)
        self.assertAlmostEqual(-3.0562684035569383, mus, 15)
        rw.reweight(lnpi, delta_mu=mus - mu)
        sat = rw.saturation_properties(volume=8**3, lnpi=lnpi, mus=mus)
        #print(sat)
        self.assertEqual(166, sat["phase_boundary"])
        self.assertAlmostEqual(0.07722558274659058, sat['vapor']['pressure'], 15)
        self.assertAlmostEqual(sat['vapor']['pressure'], sat['liquid']['pressure'], 8)
        self.assertAlmostEqual(0.10035103189192143, sat['vapor']['density'], 15)
        self.assertAlmostEqual(0.5631867734774637, sat['liquid']['density'], 15)

        #import matplotlib.pyplot as plt
        #plt.plot(lnpi)
        #plt.show()

if __name__ == "__main__":
    unittest.main()
