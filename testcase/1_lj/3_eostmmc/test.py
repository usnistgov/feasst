"""
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
"""

import argparse
import unittest
import feasst

parser = argparse.ArgumentParser()
parser.add_argument("--openMP", help="use openMP parallelization",
                    dest='openMP', action='store_true')
parser.set_defaults(openMP=False)
parser.add_argument("--boxl", "-l", help="box length", default=8, type=float)
parser.add_argument("--rCut", "-r", help="cutoff distance",
                    default=3., type=float)
parser.add_argument("--temp", "-t", help="temperature", default=1.5, type=float)
parser.add_argument("--lnz", "-z", help="activity",
                    default=-1.568214, type=float)
parser.add_argument("--nfreq", "-f", help="number of trials per print",
                    default=int(1e4), type=int)
parser.add_argument("--ncfreq", "-c", help="number of trials per print",
                    default=int(1e6), type=int)
parser.add_argument("--nMolMax", "-x", help="maximum number of mols",
                    default=5, type=int)
parser.add_argument("--molName", "-m", help="molecule file name",
                    default="data.lj", type=str)
args = parser.parse_args()
print("#", args)

def compareEnergyAndMacro(criteria, iMacro, testobject, peAv, peStd, lnPIav, lnPIstd):
    tol = 2.576*(criteria.pe(iMacro).blockStdev() + peStd)
    testobject.assertAlmostEqual(peAv, criteria.pe(iMacro).average(), delta=tol)

class TestLJ_SRSW_EOSTMMC(unittest.TestCase):
    def test(self):
        feasst.ranInitByDate()      # initialize random number generator
        space = feasst.Space(3, 0)  # initialize simulation domain
        space.lset(args.boxl)
        addMolType = space.install_dir() + "/forcefield/" + args.molName
        space.addMolInit(addMolType)

        # initialize pair-wise interactions
        pair = feasst.PairLJ(space, args.rCut)
        pair.initEnergy()

        # acceptance criteria
        nMolMin = 0
        import math
        criteria = feasst.CriteriaWLTMMC(1./args.temp, math.exp(args.lnz),
                                         "nmol", nMolMin - 0.5, args.nMolMax \
                                         + 0.5, args.nMolMax - nMolMin + 1)
        criteria.collectInit()
        criteria.tmmcInit()

        # initialize monte carlo
        mc = feasst.WLTMMC(space, pair, criteria)
        mc.weight = 3./4.
        feasst.transformTrial(mc, "translate")
        mc.weight = 1./8.
        feasst.deleteTrial(mc)
        mc.weight = 1./8.
        feasst.addTrial(mc, addMolType)

        # output log, lnpi and movie
        mc.initLog("log", args.nfreq)
        mc.initColMat("colMat", args.ncfreq)
        mc.setNFreqCheckE(args.ncfreq, 1e-8)
        mc.setNFreqTune(args.nfreq)
        mc.initMovie("movie", args.nfreq)
        #mc.initXTC("movie", args.nfreq)
        mc.initRestart("tmp/rst", args.ncfreq)

        # production tmmc simulation
        if args.openMP:
            mc.initWindows(2.,  # exponent that determines size of windows
                           0)   # extra macrostate overlap between processors
        mc.runNumSweeps(20,  # number of "sweeps"
                        -1)  # maximum number of trials. Infinite if "-1".

        # test against SRSW values
        self.assertAlmostEqual(criteria.pe(0).average(), 0, 13)
        compareEnergyAndMacro(criteria, 1, self,
                              -0.0006057402333333332,
                              6.709197666659334e-10,
                              -270.0061768+274.6763737666667,
                              0.037092307087640365)

        compareEnergyAndMacro(criteria, 2, self,
                              -0.030574223333333334,
                              9.649146611661053e-06,
                              -266.0191155333334+274.6763737666667,
                              0.03696447428346385)

        compareEnergyAndMacro(criteria, 3, self,
                              -0.089928316,
                              0.0001387472078025413,
                              -262.4277240666667+274.6763737666667,
                              0.037746391500313385)

        compareEnergyAndMacro(criteria, 4, self,
                              -0.1784570533333333,
                              3.3152449884326804e-05,
                              -259.11444086666665+274.6763737666667,
                              0.03809721387875822)

        compareEnergyAndMacro(criteria, 5, self,
                              -0.29619201333333334,
                              1.3487910636322294e-05,
                              -256.0144809+274.6763737666667,
                              0.03845757460933292)

if __name__ == "__main__":
    unittest.main()
