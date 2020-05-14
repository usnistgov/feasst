/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include <gtest/gtest.h>
#include "pair_patch_kf_multi.h"

TEST(PairPatchKFMulti, AnalyticalOnePatch) {
  auto space = feasst::makeSpace({{"dimen", "2"}, {"boxLength", "12."}});
  auto pair = feasst::makePairPatchKFMulti(space.get(), {{"rCut", "2"}, {"patchAngle", "90"}});
  pair->initData("../forcefield/data.onePatch");
  pair->rCutijset(0, 0, 1.);
  pair->rCutijset(0, 1, 0.);
  pair->rCutijset(1, 1, 2.);
  vector<double> xAdd(space->dimen());
  pair->addMol(xAdd);
  pair->initEnergy();

  // test hard sphere interaction
  xAdd[0] = 0.999;
  pair->addMol(xAdd);
  pair->initEnergy();
  EXPECT_GT(1e299, pair->peTot());
  xAdd[0] = 1.0000001 - space->x(2, 0);
  space->transMol(1, xAdd);
  pair->initEnergy();
  EXPECT_EQ(0, pair->peTot());

  // rotate the particle on the origin
  EXPECT_EQ(2, space->qMol().size());
  EXPECT_EQ(0., space->qMol()[0]);
  space->qMolAlt(0, 0, feasst::PI/2.+0.0000001);
  space->quat2pos(0);

  // test just beyond rCut for patch
  xAdd[0] = -2.000000001 - space->x(2, 0);
  space->transMol(1, xAdd);
  pair->initEnergy();
  EXPECT_EQ(0, pair->peTot());

  // test just within rCut for patch
  xAdd[0] = -1.999999 - space->x(2, 0);
  space->transMol(1, xAdd);
  pair->initEnergy();
  EXPECT_EQ(-1, pair->peTot());

  // rotate just next to but not quite on the patch angle
  space->qMolAlt(0, 0, feasst::PI/2.-0.0000001);
  space->quat2pos(0);

  pair->initEnergy();
  EXPECT_EQ(0, pair->peTot());
}

TEST(PairPatchKFMulti, AnalyticalOverlapPatch) {
  auto space = feasst::makeSpace(
   {{"dimen", "2"},
    {"boxLength", "12."}});
  auto pair = feasst::makePairPatchKFMulti(space.get(),
   {{"rCut", "2"},
    {"patchAngle", "60"}});
  // Note: rCut is the patch well distance from center
  pair->initData("../forcefield/data.overlapPatch");
  pair->initIJ();

  // modify the second patch to have a smaller patch angle
  pair->initPatchAngleInDegrees(30, 2);

  // turn off the interaction between the inner and outer patches
  pair->epsijset(1, 2, 0.);

  // add two particles
  vector<double> xAdd(space->dimen());
  pair->addMol(xAdd);
  xAdd[0] = 1.5;
  pair->addMol(xAdd);

  // rotate the second particle by 180
  vector<double> angles = {0, 0};
  EXPECT_EQ(space->qMol(), angles);
  space->qMolAlt(1, 0, feasst::PI);
  space->quat2pos(1);
  pair->initEnergy();
  EXPECT_EQ(-2, pair->peTot());

  // rotate the second particle by 135
  space->qMolAlt(1, 0, feasst::PI*0.75);
  space->quat2pos(1);
  pair->initEnergy();
  EXPECT_EQ(-1, pair->peTot());

  // rotate the second particle by 120-epsilon
  space->qMolAlt(1, 0, feasst::PI*2./3.-0.0001);
  space->quat2pos(1);
  pair->initEnergy();
  EXPECT_EQ(0, pair->peTot());
}

TEST(PairPatchKFMulti, AnalyticalOverlappingMultiSite) {
  auto space = feasst::makeSpace(
   {{"dimen", "2"},
    {"boxLength", "20."}});
  auto pair = feasst::makePairPatchKFMulti(space.get(),
   {{"rCut", "2"},
    {"patchAngle", "60"}});
  // Note: rCut is the patch well distance from center
  std::stringstream ss;
  ss << space->install_dir() << "/forcefield/data.cg4_3patch_overlap";
  pair->initData(ss.str().c_str());
  ss.str("");
  ss << space->install_dir() << "/forcefield/data.cg2_2patch_overlap_linear";
  pair->initData(ss.str().c_str());
  pair->initIJ();

  // modify the "second" duplicate patches to have a smaller patch angle
  const double innerAngle = 5; // degrees
  pair->initPatchAngleInDegrees(innerAngle, 3);
  pair->initPatchAngleInDegrees(innerAngle, 7);

  // turn off the interaction between the inner and outer patches
  pair->epsijset(2, 3, 0.);
  pair->epsijset(6, 7, 0.);
  pair->epsijset(2, 7, 0.);
  pair->epsijset(3, 6, 0.);

  // add two particles
  vector<double> xAdd(space->dimen());
  pair->addMol(xAdd);
  xAdd[1] = 3.;
  pair->addMol(xAdd);

  // also test clusters
  space->addTypeForCluster(0);
  space->addTypeForCluster(4);

  // rotate the second particle by 180
  vector<double> angles = {0, 0};
  EXPECT_EQ(space->qMol(), angles);
  space->qMolAlt(1, 0, feasst::PI);
  space->quat2pos(1);

  // test hard sphere interaction
  xAdd[1] = 2.9999999 - space->x(10, 1);
  space->transMol(1, xAdd);
  pair->initEnergy();
  EXPECT_GT(1e299, pair->peTot());
  pair->updateClusters(0.);
  EXPECT_EQ(1, space->nClusters());
  // Note: even with HS interaction, patches are still within well

  // test just outside of hs
  xAdd[1] = 3.0000001 - space->x(10, 1);
  space->transMol(1, xAdd);
  pair->initEnergy();
  EXPECT_EQ(-2, pair->peTot());
  pair->updateClusters(0.);
  EXPECT_EQ(1, space->nClusters());

  // test just inside of rCut
  xAdd[1] = 3.999999999 - space->x(10, 1);
  space->transMol(1, xAdd);
  pair->initEnergy();
  EXPECT_EQ(-2, pair->peTot());
  pair->updateClusters(0.);
  EXPECT_EQ(1, space->nClusters());

  // test turning off like-molecule interactions
  pair->epsijset(2, 2, 0.);
  pair->epsijset(3, 3, 0.);
  pair->initEnergy();
  EXPECT_EQ(0, pair->peTot());
  pair->updateClusters(0.);
  EXPECT_EQ(2, space->nClusters());

  // test just outside of rCut, and turning back on like-molecule interactions
  pair->epsijset(2, 2, 1.);
  pair->epsijset(3, 3, 1.);
  xAdd[1] = 4.000000001 - space->x(10, 1);
  space->transMol(1, xAdd);
  pair->initEnergy();
  EXPECT_EQ(0, pair->peTot());
  pair->updateClusters(0.);
  EXPECT_EQ(2, space->nClusters());

  // turn off like-molecule interactions
  pair->epsijset(2, 2, 0.);
  pair->epsijset(3, 3, 0.);
  pair->epsijset(6, 6, 0.);
  pair->epsijset(7, 7, 0.);

  // delete second particle and add cg2
  vector<int> mpart = space->imol2mpart(1);
  pair->delPart(mpart);
  space->delPart(mpart);
  xAdd[1] = 3.4;
  pair->addMol(xAdd, ss.str());
  pair->initEnergy();
  EXPECT_EQ(-2, pair->peTot());
  pair->updateClusters(0.);
  EXPECT_EQ(1, space->nClusters());

  // delete first particle and add cg2
  mpart = space->imol2mpart(0);
  pair->delPart(mpart);
  space->delPart(mpart);
  xAdd[1] = 1;
  pair->addMol(xAdd, ss.str());
  pair->initEnergy();
  EXPECT_EQ(0, pair->peTot());
  pair->updateClusters(0.);
  EXPECT_EQ(2, space->nClusters());

  // turn on like-molecule interactions and test again
  pair->epsijset(2, 2, 1.);
  pair->epsijset(3, 3, 1.);
  pair->epsijset(6, 6, 1.);
  pair->epsijset(7, 7, 1.);
  pair->initEnergy();
  EXPECT_EQ(-2, pair->peTot());
  pair->updateClusters(0.);
  EXPECT_EQ(1, space->nClusters());
  pair->printXYZ("hi", 1);

  // Analyze resulting statistics of all clusters
  EXPECT_EQ(1, space->clusterSizes().size());
  EXPECT_EQ(2, space->clusterSizes()[0]);
  // Average size of clusters, average over each configuration
  EXPECT_EQ(space->clusterSizeAccVec().vec(2).average(), 1.625);
  // Number of clusters, averaged over each configuration
  EXPECT_EQ(space->clusterNumAccVec().vec(2).average(), 1.375);
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
  EXPECT_EQ(-2, pair->peTot());
  pair->updateClusters(0.);
  EXPECT_EQ(1, space->nClusters());
  space->randRotateMulti(space->listAtoms(), -1, pair->sig());
  pair->printXYZ("movie", 0);
  space->randDispNoWrap(space->listAtoms(), -1);
  pair->printXYZ("movie", 0);
  pair->initEnergy();
  EXPECT_EQ(-2, pair->peTot());
  pair->updateClusters(0.);
  EXPECT_EQ(1, space->nClusters());
}

