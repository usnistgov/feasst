
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
import math
import feasst

class TestPairJaglaFermi_NVTMC(unittest.TestCase):
  def test(self):
    feasst.ranInitByDate()

    nMolMax=200
    space = feasst.makeSpace(feasst.args(
     {"dimen" : "2",
      "boxLength" : "20.0"}))

    #pair = feasst.PairTabular1D(space)
    pair = feasst.PairFermiJagla(space, feasst.args({"rCut" : "1.08"}))
    ssA = space.install_dir() + "/forcefield/data.fermi_jagla_a"
    pair.initData(ssA)
    ssB = space.install_dir() + "/forcefield/data.fermi_jagla_b"
    pair.initData(ssB)

    # set the cut-off
    pair.equateRcutForAllTypes()

    # set the pairwise B_0 parameter
    pair.epsijset(0, 0, 0.) # i, i, B0_ii
    pair.epsijset(0, 1, 1.) # i, j, B0_ij
    pair.epsijset(1, 1, 0.) # j, j, B0_jj

    # turn on the cell list
    space.updateCells(pair.rCut())

    # intialize the energy of the configuration
    pair.initEnergy()

    # acceptance criteria
    criteria = feasst.makeCriteriaMetropolis(feasst.args(
      {"beta" : str(10.0), # beta=1/kT
       "activ" : str(math.exp(-1))}))
    criteria.addActivity(math.exp(-1.))
    mc = feasst.MC(pair, criteria)

    # single particle translation
    mc.weight = 5.0
    feasst.addTrialTransform(mc, feasst.args(
      {"transType" : "translate",
       "maxMoveParam" : "0.1"}))

    # GCA
    mc.weight = 1./float(nMolMax)
    feasst.gcaTrial(mc)
    
    # rigid cluster moves
    space.addTypeForCluster(0)
    space.addTypeForCluster(1)
    #space.preMicellarAgg(5) # this determines the "free monomer" (freemon) stat
    mc.weight = 1./5./float(nMolMax)
    clusterCut = pair.rCut()
    feasst.clusterTrial(mc, "clustertrans", clusterCut)
    feasst.clusterTrial(mc, "clusterrotate", clusterCut)

    # position swap (expect low acceptance for strong unlike-particle ixn)
    mc.weight = 0.1
    feasst.xswapTrial(mc)

    # initialize AVB for A-B interactions only
    # NOTE: the data files also specify the AVB moves.
    # data.fermi_jagla_a specifically targets the avb of the second particle
    # type (*_b) and vice verse
    pair.neighTypeSet(0)
    pair.neighTypeSet(1)
    insDelFlag = 0  # 0: no gcinsdel, 1: gcinsdel
    dualCut = 2     # keep track of neighbors but don't use dual cut
    feasst.initConfigBias(mc, ssA, insDelFlag, dualCut)
    feasst.initConfigBias(mc, ssB, insDelFlag, dualCut)

    # initialize the number of particles
    assert(nMolMax%2==0) # assumes equimolar
    mc.nMolSeek(int(nMolMax/2), ssA)
    mc.nMolSeek(nMolMax, ssB)

    nPrint = int(1e6)
    mc.initLog("log", nPrint)
    mc.initMovie("movie", nPrint)
    mc.initRestart("tmp/rst", nPrint)
    mc.setNFreqCheckE(nPrint, 1e-6);
    mc.setNFreqTune(nPrint)
    mc.runNumTrials(int(1e7))

if __name__ == "__main__":
    unittest.main()

