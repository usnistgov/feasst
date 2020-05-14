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

class TestPairPatchKFMulti_overlap_NVTMC(unittest.TestCase):
  def test(self):
    feasst.ranInitByDate()
    space = feasst.makeSpace(feasst.args(
      {"dimen" : "2",
       "boxLength" : "40."}))
    pair = feasst.makePairPatchKFMulti(space, feasst.args(
      {"rCut" : "1.5",  # Note: rCut is the patch well center-center distance
       "patchAngle" : "31"}))
    # Note: data file sets the bond length, epsilon (depth) and sigma
    # Pair Coeffs are in order [type] [eps] [sig]
    # Patch unit vectors have 0 sigma
    # Center of mass rotation centers of 0 eps (only in cg2_2patch)
    ssA = space.install_dir() + "/forcefield/data.cg4_3patch_overlap"
    pair.initData(ssA)
    ssB = space.install_dir() + "/forcefield/data.cg2_2patch_overlap_linear"
    pair.initData(ssB)

    # If you change the sigmas, rCut is not set dependent on sigma, so this
    # initIJ may need to be changed if you want rcut-sigma=constant
    pair.initIJ() # Note: initCuts depreciated, because assigns patch angles

    """ Note that, for a pair of models, the site types are as follows

     first model:  0 - center
                   1 - hard sphere
                   2 - outer patch
                   3 - inner patch

     second model: 4 - center
                   5 - hard sphere
                   6 - outer patch
                   7 - inner patch
    """

    # modify the "second" duplicate patches to have a smaller patch angle
    # Note: this must come after "initIJ" or it will be overwritten
    innerAngle = 5 # degrees
    pair.initPatchAngleInDegrees(innerAngle, 3)
    pair.initPatchAngleInDegrees(innerAngle, 7)

    # turn off the patchy interactions between like molecules
    pair.epsijset(2, 2, 0.)
    pair.epsijset(3, 3, 0.)
    pair.epsijset(6, 6, 0.)
    pair.epsijset(7, 7, 0.)

    # turn off the interaction between the inner and outer patches
    pair.epsijset(2, 3, 0.)
    pair.epsijset(6, 7, 0.)
    pair.epsijset(2, 7, 0.)
    pair.epsijset(3, 6, 0.)

    # set the strength of the "inner" patch interaction
    pair.epsijset(3, 7, 1.)

    # set the strength of the "outer" patch interaction
    pair.epsijset(2, 6, 0.2)

    # set the minimal strength of interaction to be considered a cluster
    pair.clusterTolerance = 0.5

    # turn on cell list
    space.updateCells(pair.rCut())

    # Acceptance criteria
    criteria = feasst.makeCriteriaMetropolis(feasst.args(
      {"beta" : str(1./0.05), # beta=1/kT
       "activ" : str(math.exp(-1))}))
    criteria.addActivity(math.exp(-1.))
    mc = feasst.MC(pair, criteria)
    feasst.addTrialTransform(mc, feasst.args(
      {"transType" : "translate",
       "maxMoveParam" : "0.1"}))
    feasst.addTrialTransform(mc, feasst.args(
      {"transType" : "rotate",
       "maxMoveParam" : "0.1"}))

    # begin cluster moves
    nMolMax = 48
    mc.weight = 1./float(nMolMax)
    feasst.gcaTrial(mc)
    space.addTypeForCluster(0)
    space.addTypeForCluster(4)

    # Set the size of "premicellar aggregates"
    # This will print the concentration of "clusters" which are of size <= this value
    space.preMicellarAgg(5)

    mc.weight = 1./5./float(nMolMax)
    feasst.clusterTrial(mc, "clustertrans")
    feasst.clusterTrial(mc, "clusterrotate")

    self.assertEqual(nMolMax % 2, 0)
    mc.nMolSeek(int(round(nMolMax/2)), ssA)
    mc.nMolSeek(nMolMax, ssB)

    # Note that the movie file outputs "patches" as inscribed spheres.
    # From PairPatchKFMulti::printxyz line 179:
    # make a patch by inscribing sphere of radius r2 inside bead, radius r1=sig
    #  distance between center of sphere and bead is a
    #  for given patch angle, use law of cosines to derive r2
    #  constraint r1 + eps = r2 + a, where eps is small, to avoid clipping
    nPrint = int(1e4)
    mc.initLog("log", nPrint)
    mc.initMovie("movie", nPrint)
    mc.initRestart("tmp/rst", nPrint)
    mc.setNFreqCheckE(nPrint, 1e-6)
    mc.setNFreqTune(nPrint)

    # reset cluster statistics after equilibration
    space.clusterReset()

    mc.runNumTrials(int(1e6))

if __name__ == "__main__":
    unittest.main()
