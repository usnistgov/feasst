#include <gtest/gtest.h>
#include "pair_patch_kf.h"
#include "mc.h"
#include "trial_transform.h"

using namespace feasst;

TEST(PairPatchKF, patchKFAnalytical1) {
  Space s(3,0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(10,dim);
  s.readXYZBulk(2, "onePatch", "../unittest/patch/onePatch6.xyz");
  PairPatchKF p(&s, 3, 90);
  EXPECT_NEAR(0, p.cpa(), 1e-15);
  p.initEnergy();
  EXPECT_NEAR(0, p.peTot(), 1e-15);
}

TEST(PairPatchKF, patchKFAnalytical2) {
  Space s(3,0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(10,dim);
  s.readXYZBulk(2, "onePatch", "../unittest/patch/onePatch5.xyz");
  PairPatchKF p(&s, 3, 90);
  EXPECT_NEAR(0, p.cpa(), 1e-15);
  p.initEnergy();
  const double petot = p.peTot();
  EXPECT_NEAR(-3, petot, 1e-15);
  EXPECT_NEAR(petot, p.allPartEnerForce(0), DTOL);
  EXPECT_NEAR(petot, p.allPartEnerForce(1), DTOL);
  EXPECT_EQ(1, p.checkEnergy(1e-18, 1));

  s.addTypeForCluster(0);
  p.updateClusters(0);
  EXPECT_EQ(3, s.nClusters());

}

TEST(PairPatchKF, patchKFcellList) {
  Space s(3,0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(10,dim);
  s.readXYZBulk(2, "onePatch", "../unittest/patch/onePatch50.xyz");
  PairPatchKF p(&s, 1.5, 90);
  p.initEnergy();
  const double petot = p.peTot();
  EXPECT_NEAR(-62, petot, 1e-15);

  // see if cell list gives same result
  s.updateCells(p.rCut(), p.rCut());
  EXPECT_EQ(1, p.checkEnergy(1e-18, 1));
}

TEST(PairPatchKF, patchKFmirrorAnalytical) {
  Space s(3,0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(10,dim);
  s.readXYZBulk(2, "onePatch", "../unittest/patch/onePatch5.xyz");
  PairPatchKF p(&s, 3, 90);
  p.mirrorPatch(1);
  p.printxyz("tmp/onePatch5vis.xyz", 1);
  EXPECT_NEAR(0, p.cpa(), 1e-15);
  p.initEnergy();
  EXPECT_NEAR(-7, p.peTot(), 1e-15);

  const double chi = 0.4;
  const double theta = acos(1-chi)*180./PI;
  EXPECT_NEAR(53.1301023541560000, theta, 1e-13);
  PairPatchKF p2(&s, 3, theta);
  p2.mirrorPatch(1);
  p2.printxyz("tmp/onePatch5vis.xyz", 0);
  EXPECT_NEAR(1-chi, p2.cpa(), 1e-15);
  p2.initEnergy();
  const double petot = p2.peTot();
  EXPECT_NEAR(-2, petot, 1e-15);

  EXPECT_EQ(1, p.checkEnergy(1e-18, 1));
}

TEST(MC, PairPatchKF) {
  feasst::Space space;
  space.lset(8);
  space.addMolInit("../forcefield/data.onePatch");
  feasst::PairPatchKF pair(&space, 2., 90);
  pair.initData("../forcefield/data.onePatch");
  pair.initEnergy();
  feasst::CriteriaMetropolis criteria(1., exp(-1));
  feasst::MC mc(&space, &pair, &criteria);
  feasst::transformTrial(&mc, "translate", 1);
  feasst::transformTrial(&mc, "rotate", 1);
//  space.addTypeForCluster(0);
//  feasst::clusterTrial(&mc, "clustertrans", -1, 1.);
//  feasst::clusterTrial(&mc, "clusterrotate", -1, 1.);
  // feasst::gcaTrial(&mc);
  mc.nMolSeek(12, "../forcefield/data.onePatch");
  pair.initEnergy();
  const int nfreq = 1e0;
  mc.initLog("tmp/patchlog", nfreq);
  mc.initMovie("tmp/patchmovie", nfreq);
  mc.setNFreqCheckE(nfreq, feasst::DTOL);
  mc.runNumTrials(1e2);
}

