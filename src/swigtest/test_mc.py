"""
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
"""

"""@docstring
@brief Tests Monte Carlo
@author Nathan A. Mahynski
@date 5/16/2017
@filename test_mc.py
"""

import sys, os
import unittest
import numpy as np

class TestMC(unittest.TestCase):
    def testMC(self):
        s=feasst.Space(3,0)
        s.lset(10)
        s.addMolInit("../forcefield/data.lj")
        p=feasst.PairLJMulti(s, 3)
        p.initData("../forcefield/data.lj")
        p.cutShift(1)
        c=feasst.CriteriaWLTMMC(1, 0.1, "nmol",-0.5,5.5,6)
        mc=feasst.WLTMMC(s, p, c)
        feasst.transformTrial(mc, "translate")
        feasst.addTrial(mc, "../forcefield/data.lj")
        feasst.deleteTrial(mc)
        mc.runNumTrials(100)
        self.assertFalse(mc.peAccumulator().sumDble() == 0)
        mc.peAccumulator().reset()
        mc.zeroStat()
        self.assertTrue(mc.peAccumulator().sumDble() == 0)

if __name__ == "__main__":
    def find_all(name, path):
        result = []
        for root, dirs, files in os.walk(path):
            if name in files:
                result.append(os.path.join(root, name))
        return result

    feasstname = find_all("_feasst.so", os.path.dirname(os.path.realpath(__file__))+"/../../")
    if len(feasstname) > 1:
        print ('Found multiple FEASST sources which may cause an error: ', feasstname)
        print ("using first one")
    feasstdir = feasstname[0].split('/_feasst.so')[0]
    sys.path.append(feasstdir)

    import feasst
    unittest.main()
