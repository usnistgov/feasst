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

class TestLJ_SLITPORE_EOSTMMC(unittest.TestCase):
  def test(self):
    feasst.ranInitByDate()
    space = feasst.makeSpace(feasst.args(
     {"dimen" : "3",
      "boxLength" : "10"}))
    confine_dim = 2  # apply confinement in the z-dimension
    upper = space.boxLength(0)/2.
    lower = -upper
    space.initBoxLength(space.boxLength(1) + 0, confine_dim)     # stretch PBC box length for slab
    space.initPBC(confine_dim, 0) # disable periodicity(0) in confined dimension

    # Initialize LJ interactions
    pairLJ = feasst.makePairLJ(space, feasst.args(
     {"rCut" : "3.",
      "cutType" : "lrc",
      "molTypeInForcefield" : "data.lj"}))
    pairLJ.initEnergy()

    # Initialize Wall
    pairWall1 = feasst.makePairFieldSlitLJ(space)
    pairWall1.initSlit(confine_dim, upper, lower)
    pairWall2 = pairWall1.clone(space)
    pairWall1.initAlpha(9)
    pairWall1.initEps(0, 2./15.)
    pairWall1.initSig(0, 1.25)
    pairWall2.initAlpha(3)
    pairWall2.initEps(0, -1.)
    pairWall2.initSig(0, 1.2)
    pairWall2.initDelta(0.25)

    # Create PairHybrid object which encompasses both LJ and wall interaction
    pair = feasst.makePairHybrid(space)
    pair.addPair(pairWall1)  # add cheaper pairs first for efficiency
    pair.addPair(pairWall2)
    pair.addPair(pairLJ)
    pair.initEnergy()

    # Initialize acceptance criteria and flat-histogram
    import math
    activ = math.exp(-11)
    temp = 1.5
    nMolMax = 100
    criteria = feasst.makeCriteriaWLTMMC(feasst.args(
     {"beta" : str(1./temp),
      "activ" : str(activ),
      "mType" : "nmol",
      "nMax" : str(nMolMax)}))
    criteria.collectInit()  # initialize updating of collection matrix
    criteria.tmmcInit()     # initialize use of collection matrx (e.g., TMMC)

    # Initialize WLTMMC and trial moves
    mc = feasst.WLTMMC(space, pair, criteria)

    mc.weight = 0.4
    feasst.transformTrial(mc, "translate")
    feasst.transformTrial(mc, "rotate")

    tdel = feasst.TrialDelete(space.addMolListType(0))
    # tdel.numFirstBeads(10)
    mc.weight = 0.1
    mc.initTrial(tdel)

    tadd = feasst.TrialAdd(space.addMolListType(0))
    # tadd.numFirstBeads(10)
    mc.weight = 0.1
    mc.initTrial(tadd)

# Note, this can be turned on or off and shoudln't affect the result
# Perhaps turning it on will give a slight efficiency boost
#    # Note that confine can't be called until the space object pointer is
#    # initialized. This means it must follow "initTrial" in this example.
#    tdel.confine(upper, lower, confine_dim)
#    tadd.confine(upper, lower, confine_dim)

    # Initialize outputs and checks
    nPrint = int(1e4)
    mc.initLog("log", nPrint)
    mc.initMovie("movie", nPrint)
    # mc.initXTC("traj", nPrint)
    mc.initColMat("colMat", nPrint)
    mc.initRestart("tmp/rst", nPrint)
    mc.setNFreqCheckE(10*nPrint, 1e-4)
    mc.setNFreqTune(nPrint)

    # Initialize parallel windows by dividing order parameter, nmol
    #mc.initWindows(1)
    #mc.runNumSweeps(10, -1)

    # Serial run
    mc.initProduction()
    mc.runNumTrials(int(1e7))

    # print some statistics
    print("Average potential")
    for ipair in range(pair.nPairs()):
      ppair = pair.getPair(ipair)
      print(" pair" + str(ipair), ppair.accumulator().average(), ppair.accumulator().std())

if __name__ == "__main__":
    unittest.main()
