#include <gtest/gtest.h>
#include "pair_patch_kf.h"

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

  EXPECT_EQ(1, p.checkEnergy(1e-18, 1));
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

