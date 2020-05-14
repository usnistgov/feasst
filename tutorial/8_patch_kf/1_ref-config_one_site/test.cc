/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "feasst.h"

int main() {  // PatchKF, REFCONF_REINHARDT
  auto space = feasst::makeSpace(
   {{"dimen", "2"},
    {"boxLength", "20."}});
  auto pair = feasst::makePairPatchKFMulti(space.get(),
   {{"rCut", "2"},
    {"patchAngle", "30"}});
  // Note: rCut is the patch well distance from center
  pair->patchType = 0;
  std::stringstream ss;
  ss << space->install_dir() << "/forcefield/data.onePatch";
  pair->initData(ss.str().c_str());
  pair->initIJ();
  vector<double> xAdd(space->dimen());
  pair->addMol(xAdd);
  xAdd[0] = 0.;
  pair->addMol(xAdd);

  // also test clusters
  space->addTypeForCluster(0);

  // test hard sphere interaction
  xAdd[0] = 0.99 - space->x(2, 0);
  space->transMol(1, xAdd);
  pair->initEnergy();
  ASSERT(1e297 < pair->peTot(), "HS0");
  pair->updateClusters(0.);
  pair->printXYZ("C1", 1);
  ASSERT(2 == space->nClusters(), "C0");
  // Note: even with HS interaction, patches are still within well

  // test no interaction beyond hard sphere
  xAdd[0] = 1.000000001 - space->x(2, 0);
  space->transMol(1, xAdd);
  pair->initEnergy();
  ASSERT(pair->peTot() < feasst::DTOL, "HS0");
  pair->updateClusters(0.);
  pair->printXYZ("C1", 1);
  ASSERT(2 == space->nClusters(), "C0");
  // Note: even with HS interaction, patches are still within well

  // rotate the second particle by 180
  vector<double> angles = {0, 0};
  ASSERT(space->qMol() == angles, "Initial angles expected [0, 0]");
  space->qMolAlt(1, 0, feasst::PI);
  space->quat2pos(1);

  // test hard sphere interaction
  xAdd[0] = 0.99 - space->x(2, 0);
  space->transMol(1, xAdd);
  pair->initEnergy();
  ASSERT(1e297 < pair->peTot(), "HS1");
  pair->updateClusters(0.);
  pair->printXYZ("C1", 1);
  ASSERT(1 == space->nClusters(), "C1");
  // Note: even with HS interaction, patches are still within well

  // test just outside of hs
  xAdd[0] = 1.0000001 - space->x(2, 0);
  space->transMol(1, xAdd);
  pair->initEnergy();
  ASSERT(-1 == pair->peTot(), "SQW1");
  pair->updateClusters(0.);
  ASSERT(1 == space->nClusters(), "CSQ1");

  // test just outside of rCut
  xAdd[0] = 2.0000000001 - space->x(2, 0);
  space->transMol(1, xAdd);
  pair->initEnergy();
  ASSERT(pair->peTot() < feasst::DTOL, "SQW3");
  pair->updateClusters(0.);
  ASSERT(2 == space->nClusters(), "C3");

  // test just inside of rCut
  xAdd[0] = 1.999999999 - space->x(2, 0);
  space->transMol(1, xAdd);
  pair->initEnergy();
  ASSERT(-1 == pair->peTot(), "SQW2");
  pair->updateClusters(0.);
  ASSERT(1 == space->nClusters(), "C2");

  // Analyze resulting statistics of all clusters
  ASSERT(1 == space->clusterSizes().size(), "clustSizes");
  ASSERT(2 == space->clusterSizes()[0], "clus");
  ASSERT(space->clusterSizeAccVec().vec(2).average() == 1.5,
    "Average size of clusters, average over each configuration");
  ASSERT(space->clusterNumAccVec().vec(2).average() == 1.5,
    "Number of clusters, averaged over each configuration");
  // Note the cluster size probability distribution is defined as:
  // -  pick a random cluster. What is the probility it is of size "x"?
  space->printClusterStat("clus");

  // Displace particles such that they span the periodic boundary,
  // then attempt a rigid cluster rotation and translation
  pair->printXYZ("movie", 1);
  xAdd[1] = -space->boxLength(0)/2.;
  space->transMol(0, xAdd);
  xAdd[1] = +space->boxLength(0)/2.;
  space->transMol(1, xAdd);
  pair->printXYZ("movie", 0);
  pair->initEnergy();
  ASSERT(-1 == pair->peTot(), "PBC1");
  pair->updateClusters(0.);
  ASSERT(1 == space->nClusters(), "C7");
  space->randRotateMulti(space->listAtoms(), -1, pair->sig());
  pair->printXYZ("movie", 0);
  space->randDispNoWrap(space->listAtoms(), -1);
  pair->printXYZ("movie", 0);
  pair->initEnergy();
  ASSERT(-1 == pair->peTot(), "PBC2");
  pair->updateClusters(0.);
  ASSERT(1 == space->nClusters(), "C8");
}

