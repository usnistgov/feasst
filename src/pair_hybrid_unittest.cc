
#include <gtest/gtest.h>
#include "pair_lj.h"
#include "pair_lj_multi.h"
#include "pair_lj_coul_ewald.h"
#include "pair_patch_kf.h"
#include "pair_hybrid.h"

TEST(PairHybrid, hybrid) {
  Space s(3,0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(24.8586887,dim);
  s.readXYZBulk(3, "water", "../unittest/spce/test52.xyz");
  s.addMolInit("../forcefield/data.spce");
  PairLJCoulEwald p(&s, 12.42934435);
  p.initBulkSPCE(5.6, 38.00000001);
  double petot = p.peTot();
  PairHybrid ph(&s, 12);
  ph.addPair(&p);
  EXPECT_EQ(1, ph.nPairs());
  EXPECT_EQ(petot, ph.peTot());

  Space s2(3, 0);
  for (int dim=0; dim < s2.dimen(); ++dim) s2.lset(8,dim);
  s2.readXYZBulk(1, "atom", "../unittest/lj/srsw/lj_sample_config_periodic4.xyz");
  PairLJ p2(&s2, 3);
  p2.linearShift(1);
  p2.initEnergy();
  petot += p2.peTot();
  ph.addPair(&p2);
  EXPECT_EQ(2, ph.nPairs());
  EXPECT_EQ(petot, ph.peTot());

  Space s3(3,0);
  for (int dim=0; dim < s3.dimen(); ++dim) s3.lset(10,dim);
  s3.readXYZBulk(2, "onePatch", "../unittest/patch/onePatch5.xyz");
  PairPatchKF p3(&s3, 3, 90);
  p3.initEnergy();
  petot += p3.peTot();
  ph.addPair(&p3);
  EXPECT_EQ(3, ph.nPairs());
  EXPECT_EQ(petot, ph.peTot());

}

TEST(PairHybrid, ljANDwca) {

  Space s(3, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(60,dim);
  s.addMolInit("../forcefield/data.cg3_60_43_1");
  s.addMol("../forcefield/data.cg3_60_43_1");
  s.addMol("../forcefield/data.cg3_60_43_1");
  s.addMol("../forcefield/data.cg3_60_43_1");
  s.addMol("../forcefield/data.cg3_60_43_1");
  s.addMol("../forcefield/data.cg3_60_43_1");
  PairLJMulti pLJ(&s, 3);
  pLJ.initLMPData("../forcefield/data.cg3_60_43_1");
  pLJ.linearShift(1);
  pLJ.epsijset(1, 2, 0.);
  pLJ.epsijset(2, 1, 0.);
  pLJ.epsijset(2, 2, 0.);
  pLJ.initEnergy();
  PairLJMulti pWCA(&s, pow(2, 1./6.));
  pWCA.initLMPData("../forcefield/data.cg3_60_43_1");
  pWCA.cutShift(1);
  pWCA.epsijset(1, 1, 0.);
  pWCA.initEnergy();
  PairHybrid ph(&s, s.minl()/2.);
  ph.addPair(&pLJ);
  ph.addPair(&pWCA);
  EXPECT_EQ(2, ph.nPairs());
  //EXPECT_NEAR(0.752213404790411000, ph.peTot(), 1e-16);
  EXPECT_EQ(1, pLJ.checkEnergy(1e-10, 1));
  EXPECT_EQ(1, pWCA.checkEnergy(1e-10, 1));
  EXPECT_EQ(1, ph.checkEnergy(1e-10, 1));
  EXPECT_EQ(1, pLJ.checkEnergy(1e-10, 0));

  //cout << "about to move" << endl;
  vector<int> mpart(4);
  mpart[0] = 0; mpart[1] = 1; mpart[2] = 2; mpart[3] = 3;
  ph.multiPartEner(mpart, 0);
  ph.update(mpart, 0, "store");
  s.randDisp(mpart, 2);
  s.randRotate(mpart, 2);
  ph.multiPartEner(mpart, 1);
  ph.update(mpart, 0, "update");
  //cout << "moved" << endl;
  EXPECT_EQ(1, pLJ.checkEnergy(1e-10, 0));
  EXPECT_EQ(1, pWCA.checkEnergy(1e-10, 0));
  EXPECT_EQ(1, ph.checkEnergy(1e-10, 0));

//  cout << "about to add" << endl;
//  s.addMol("../forcefield/data.cg3_60_43_1");
//  ph.addPart();
//  cout << "added" << endl;
//  EXPECT_EQ(1, pLJ.checkEnergy(1e-10, 0));
//  EXPECT_EQ(1, pWCA.checkEnergy(1e-10, 0));
//  EXPECT_EQ(1, ph.checkEnergy(1e-10, 0));

//  EXPECT_EQ(1, pWCA.checkEnergy(1e-10, 0));
//  EXPECT_EQ(1, ph.checkEnergy(1e-10, 0));
//  EXPECT_EQ(1, pLJ.checkEnergy(1e-10, 1));
//  EXPECT_EQ(1, pWCA.checkEnergy(1e-10, 1));
//  EXPECT_EQ(1, ph.checkEnergy(1e-10, 1));

  // clone
  shared_ptr<Space> s2 = s.cloneShrPtr();
  Pair* p2 = ph.clone(s2.get());
  delete p2;

  // restart
  s.writeRestart("tmp/srst");
  Space s3("tmp/srst");
  pLJ.writeRestart("tmp/pLJrst");
  PairLJMulti pLJ3(&s, "tmp/pLJrst");
  pWCA.writeRestart("tmp/pWCArst");
  PairLJMulti pWCA3(&s, "tmp/pWCArst");
  ph.writeRestart("tmp/phrst");
  PairHybrid ph3(&s, "tmp/phrst");

}

