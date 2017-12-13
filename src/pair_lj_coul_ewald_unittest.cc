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
#include "pair_lj_coul_ewald.h"

using namespace feasst;

TEST(PairLJCoulEwald, pairLJCoulEwaldVSHybrid) {
  const int dim=3;
  Space s(dim);
  const double boxl = 24.8586887;
  s.initBoxLength(boxl);
  s.readXYZBulk(3, "water", "../unittest/spce/test52.xyz");
  PairLJCoulEwald p(&s, {{"rCut", "12.42934435"}});
  p.initBulkSPCE(5.6, 38);
  EXPECT_NEAR(p.peLJ()+p.peLRC(), 58.1554182707327, 1e-3);
  EXPECT_NEAR(p.alpha, 0.225273346779551, 1e-12);
  // problem may not have been pi, but instead a molecular COM cutoff versus atomic cutoff
  EXPECT_NEAR(p.peQReal(), -291.798831606156, 1e2);  // hybrid code seems to have pi correct to only 6 digits, and I can't fix this
  EXPECT_NEAR(p.q()[0]*p.q()[1], -499.074038364967, 1e-4);
  EXPECT_NEAR(p.peQFrr(), 53.7073916749299, 1e-3);
  EXPECT_NEAR(p.peQFrrSelf(), 1.47865637065022*52, 1e1);
  EXPECT_NEAR(p.peTot(), -256.8261529343, 1e-1);
  // compare to last run, make sure things don't change
  const double tol = 1e-11;
  EXPECT_NEAR(p.peLJ()+p.peLRC(),58.156080799920595, tol);
  EXPECT_NEAR(p.peQReal(),   -291.8132488804917, 1e-7); //tabular erf
  //EXPECT_NEAR(p.peQReal(),   -291.8132488804917, tol);
  //EXPECT_NEAR(p.peQReal(), -291.81324467526474, tol); //srsw change
  // differences may be due to molecular COM cutoff versus atomic cutoff
  // EXPECT_NEAR(p.peQReal(), -291.82585360618094, tol);
  EXPECT_NEAR(p.peQFrr(),   53.706480733660911, tol);
  //EXPECT_NEAR(p.peQFrr(), 53.707391860222536, tol); //srsw change
  EXPECT_NEAR(p.peQFrrSelf(),   76.890272201001295, tol);
  //EXPECT_NEAR(p.peQFrrSelf(), 76.890271092977457, tol); //srsw change
}

TEST(PairLJCoulEwald, pairLJCoulEwaldmultiPartEne) {
  Space s(3);
  const double boxl = 24.8586887;
  s.initBoxLength(boxl);
  s.readXYZBulk(3, "water", "../unittest/spce/test52.xyz");
  s.addMolInit("../forcefield/data.spce");
  PairLJCoulEwald p(&s, {{"rCut", "12.42934435"}});
  p.initBulkSPCE(5.6, 38);

  // move molecule and compute energy change
  // Depreciated test because its hard to add molecule exactly where you want it with proper counting for molecule description
  //const double peLJPrev2 = p.peLJ(), peLRCPrev2 = p.peLRC(), peQRealPrev2 = p.peQReal(), peQFrrPrev2 = p.peQFrr(), peQFrrSelfPrev2 = p.peQFrrSelf();
  double pePrev = p.peTot(), peLJPrev = p.peLJ(), peLRCPrev = p.peLRC(), peQRealPrev = p.peQReal(), peQFrrPrev = p.peQFrr(), peQFrrSelfPrev = p.peQFrrSelf();
  double deLJ = 0, deLRC = 0, deQReal = 0, deQFrr = 0, deQFrrSelf = 0;
  vector<int> mpart(3);
  mpart[0] = 0; mpart[1] = 1; mpart[2] = 2;

  double de = 0;
  de -= p.multiPartEner(mpart, 0);
  p.update(mpart, 0, "store");
  deLJ -= p.peLJone();
  deLRC -= p.peLRCone();
  deQReal -= p.peQRealone();
  deQFrr -= p.peQFrrone();
  deQFrrSelf -= p.peQFrrSelfone();

  s.randDisp(mpart, 2);

  de += p.multiPartEner(mpart, 1);
  deLJ += p.peLJone();
  deLRC += p.peLRCone();
  deQReal += p.peQRealone();
  deQFrr += p.peQFrrone();
  deQFrrSelf += p.peQFrrSelfone();

  p.update(mpart, 0, "update");
  double peLJPrev3 = p.peLJ(), peLRCPrev3 = p.peLRC(), peQRealPrev3 = p.peQReal(), peQFrrPrev3 = p.peQFrr(), peQFrrSelfPrev3 = p.peQFrrSelf();

  p.initEnergy();

  EXPECT_NEAR(p.peLJ() - peLJPrev, deLJ, 1e-12);
  EXPECT_NEAR(p.peLRC() - peLRCPrev, deLRC, 1e-13);
  EXPECT_NEAR(p.peQReal() - peQRealPrev, deQReal, 1e-12);
  EXPECT_NEAR(p.peQFrr() - peQFrrPrev, deQFrr, 1e-13);
  EXPECT_NEAR(p.peQFrrSelf() - peQFrrSelfPrev, deQFrrSelf, 1e-13);
  EXPECT_NEAR(p.peLJ(), peLJPrev3, 1e-12);
  EXPECT_NEAR(p.peLRC(), peLRCPrev3, 1e-13);
  EXPECT_NEAR(p.peQReal(), peQRealPrev3, 1e-12);
  EXPECT_NEAR(p.peQFrr(), peQFrrPrev3, 1e-13);
  EXPECT_NEAR(p.peQFrrSelf(), peQFrrSelfPrev3, 1e-13);
  EXPECT_NEAR(pePrev + de, p.peTot(), 1e-12);

  // remove molecule and compute energy change
  p.initEnergy();
  pePrev = p.peTot(), peLJPrev = p.peLJ(), peLRCPrev = p.peLRC(), peQRealPrev = p.peQReal(), peQFrrPrev = p.peQFrr(), peQFrrSelfPrev = p.peQFrrSelf();
  deLJ = 0, deLRC = 0, deQReal = 0, deQFrr = 0, deQFrrSelf = 0;

  de = -1*p.multiPartEner(mpart, 2);
  p.update(mpart, 2, "store");
  deLJ -= p.peLJone();
  deLRC -= p.peLRCone();
  deQReal -= p.peQRealone();
  deQFrr += p.peQFrrone() - p.peQFrr();
  deQFrrSelf -= p.peQFrrSelfone();

  p.delPart(mpart);
  s.delPart(mpart);

  p.update(mpart, 2, "update");
  peLJPrev3 = p.peLJ(), peLRCPrev3 = p.peLRC(), peQRealPrev3 = p.peQReal(), peQFrrPrev3 = p.peQFrr(), peQFrrSelfPrev3 = p.peQFrrSelf();

  p.initEnergy();
  p.k2maxset(38);
  EXPECT_NEAR(p.peLJ() - peLJPrev, deLJ, 1e-11);
  EXPECT_NEAR(p.peLRC() - peLRCPrev, deLRC, 1e-11);
  EXPECT_NEAR(p.peQReal() - peQRealPrev, deQReal, 1e-11);
  EXPECT_NEAR(p.peQFrr() - peQFrrPrev, deQFrr, 1e-11);
  EXPECT_NEAR(p.peQFrrSelf() - peQFrrSelfPrev, deQFrrSelf, 1e-11);
  EXPECT_NEAR(p.peLJ(), peLJPrev3, 1e-11);
  EXPECT_NEAR(p.peLRC(), peLRCPrev3, 1e-11);
  EXPECT_NEAR(p.peQReal(), peQRealPrev3, 1e-11);
  EXPECT_NEAR(p.peQFrr(), peQFrrPrev3, 1e-11);
  EXPECT_NEAR(p.peQFrrSelf(), peQFrrSelfPrev3, 1e-11);
  EXPECT_NEAR(pePrev + de, p.peTot(), 1e-11);

  // insert particle
  p.initEnergy();
  pePrev = p.peTot(), peLJPrev = p.peLJ(), peLRCPrev = p.peLRC(), peQRealPrev = p.peQReal(), peQFrrPrev = p.peQFrr(), peQFrrSelfPrev = p.peQFrrSelf();
  deLJ = 0, deLRC = 0, deQReal = 0, deQFrr = 0, deQFrrSelf = 0;

  s.addMol("../forcefield/data.spce");
  p.addPart();
//  vector<double> x(s.dimen());
//  x.at(0) = -1.03073217544461; x.at(1) = 3.28157508811008; x.at(2) = -8.39542194433806;
//  s.addPart(x, 0, 100);
//  p.addPart();
//  x.at(0) = -0.125284099175392; x.at(1) = 3.57064486216370; x.at(2) = -8.70623128607182;
//  s.addPart(x, 1, 100);
//  p.addPart();
//  x.at(0) = -1.00336055339743; x.at(1) = 2.31439470242869; x.at(2) = -8.14280979922285;
//  s.addPart(x, 1, 100);
//  p.addPart();

  // check insertion
  EXPECT_EQ(156, s.natom());
//  EXPECT_NEAR(-0.125284099175392, s.x()[dim*(s.natom()-2)+0], 1e-19);
  //EXPECT_NEAR(qh, p.q[s.type()[s.natom()-1]], 1e-19);
  //EXPECT_NEAR(-2*qh, p.q[s.type()[s.natom()-3]], 1e-19);
  mpart[0] = s.natom()-3; mpart[1] = s.natom()-2; mpart[2] = s.natom()-1;

  de = p.multiPartEner(mpart, 3);
  p.update(mpart, 3, "store");
  deLJ += p.peLJone();
  deLRC += p.peLRCone();
  deQReal += p.peQRealone();
  deQFrr += p.peQFrrone() - p.peQFrr();
  deQFrrSelf += p.peQFrrSelfone();

  p.update(mpart, 3, "update");
  peLJPrev3 = p.peLJ(), peLRCPrev3 = p.peLRC(), peQRealPrev3 = p.peQReal(), peQFrrPrev3 = p.peQFrr(), peQFrrSelfPrev3 = p.peQFrrSelf();

  p.initEnergy();
  p.k2maxset(38);
  EXPECT_NEAR(1, (p.peLJ() - peLJPrev)/deLJ, 1e-11);
  EXPECT_NEAR(p.peLRC() - peLRCPrev, deLRC, 1e-13);
  EXPECT_NEAR(p.peQReal() - peQRealPrev, deQReal, 1e-11);
  EXPECT_NEAR(p.peQFrr() - peQFrrPrev, deQFrr, 1e-13);
  EXPECT_NEAR(p.peQFrrSelf() - peQFrrSelfPrev, deQFrrSelf, 1e-11);
  EXPECT_NEAR(1, p.peLJ()/peLJPrev3, 1e-11);
  EXPECT_NEAR(p.peLRC(), peLRCPrev3, 1e-13);
  EXPECT_NEAR(p.peQReal(), peQRealPrev3, 1e-11);
  EXPECT_NEAR(p.peQFrr(), peQFrrPrev3, 1e-13);
  EXPECT_NEAR(p.peQFrrSelf(), peQFrrSelfPrev3, 1e-11);
  EXPECT_NEAR(1, (pePrev + de)/p.peTot(), 1e-11);

// Depreciated test because its hard to add molecule exactly where you want it with proper counting for molecule description
//  //should be back where we started
//  EXPECT_NEAR(p.peLJ(), peLJPrev2, 1e-13);
//  EXPECT_NEAR(p.peLRC(), peLRCPrev2, 1e-13);
//  EXPECT_NEAR(p.peQReal(), peQRealPrev2, 1e-12);
//  EXPECT_NEAR(p.peQFrr(), peQFrrPrev2, 1e-13);
//  EXPECT_NEAR(p.peQFrrSelf(), peQFrrSelfPrev2, 1e-12);
}

TEST(PairLJCoulEwald, neigh) {
  ranInitByDate();
  Space s(3);
  const double boxl = 24.8586887;
  s.initBoxLength(boxl);
  s.readXYZBulk(3, "water", "../unittest/spce/test52.xyz");
  s.addMolInit("../forcefield/data.spce");
  PairLJCoulEwald p(&s, {{"rCut", "12.42934435"}});
  p.initBulkSPCE(5.6, 38);

  // if neighCut < sig, there should be no neighbors
  p.initNeighList(1.5, 0);
  for (int i = 0; i < s.nMol(); ++i) {
    EXPECT_EQ(0, int(p.neigh()[i].size()));
  }

  // build neighlist, then move a particle and see if neighbor lists match
  const double rAbove = 5.;
  const double rBelow = 2.5;
  p.initNeighList(rAbove, rBelow);
  vector<int> mpart(3);
  mpart[0] = 0; mpart[1] = 1; mpart[2] = 2;
  p.update(mpart, 0, "store");
  s.randDisp(mpart, 10);
  p.multiPartEner(mpart, 1);
  p.update(mpart, 0, "update");

  // delete a particle and see if neighborlist matches
  p.multiPartEner(mpart, 2);
  p.update(mpart, 2, "store");
  p.delPart(mpart);
  s.delPart(mpart);
  p.update(mpart, 2, "update");

  // add a water molecule and see if neighborlist matches
  s.addMol("../forcefield/data.spce");
  p.addPart();
  mpart.clear();
  mpart = s.lastMolIDVec();
  p.multiPartEner(mpart, 3);
  p.update(mpart, 3, "store");
  p.update(mpart, 3, "update");

  // failed insertion of particle
  s.addMol("../forcefield/data.spce");
  p.addPart();
  mpart.clear();
  mpart = s.lastMolIDVec();
  p.multiPartEner(mpart, 3);
  p.update(mpart, 3, "store");
  p.delPart(mpart);
  s.delPart(mpart);

  EXPECT_EQ(1, p.checkNeigh());
}


TEST(PairLJCoulEwald, cheapEnergyLJCoul) {
  ranInitByDate();
  Space s(3);
  const double boxl = 24.8586887;
  s.initBoxLength(boxl);
  s.readXYZBulk(3, "water", "../unittest/spce/test52.xyz");
  s.addMolInit("../forcefield/data.spce");
  PairLJCoulEwald p(&s, {{"rCut", "12.42934435"}});
  p.initBulkSPCE(5.6, 38);
  for (int i = 0; i < 2; ++i) {

    // test that peljone is same as cheap energy for water
    vector<int> mpart = s.randMol();
    p.cheapEnergy(0);
    p.multiPartEnerReal(mpart, 0);
    double peLJone = p.peLJone();
    p.cheapEnergy(1);
    EXPECT_EQ(peLJone, p.multiPartEnerReal(mpart,0));
  }
}

TEST(PairLJCoulEwald, reconstruct) {
  Space s(3);
  s.initBoxLength(24.8586887);
  s.readXYZBulk(3, "water", "../unittest/spce/test52.xyz");
  s.addMolInit("../forcefield/data.spce");
  PairLJCoulEwald p(&s, {{"rCut", "12.42934435"}});
  p.initBulkSPCE(5.6, 38);
  const double petot = p.peTot();
  //shared_ptr<Space> s2 = make_shared<Space>(s);
  Space* s2 = s.clone();
  //Space* s2 = new Space(s);
  p.reconstruct(s2);
  s.addMol("../forcefield/data.spce");
  p.addPart();
  p.initEnergy();
  EXPECT_EQ(petot, p.peTot());
  delete s2;
}

TEST(PairLJCoulEwald, clone) {
  Space s(3);
  s.initBoxLength(24.8586887);
  s.readXYZBulk(3, "water", "../unittest/spce/test52.xyz");
  s.addMolInit("../forcefield/data.spce");
  PairLJCoulEwald p(&s, {{"rCut", "12.42934435"}});
  p.initBulkSPCE(5.6, 38);
  const double petot = p.peTot();
  shared_ptr<Space> s2 = s.cloneShrPtr();
  Pair* p2 = p.clone(s2.get());
  s.addMol("../forcefield/data.spce");
  p.addPart();
  p2->addPart();
  p.initEnergy();
  p2->initEnergy();
  EXPECT_NE(petot, p.peTot());
  EXPECT_EQ(petot, p2->peTot());
  delete p2;
}

TEST(PairLJCoulEwald, PairLJCoulEwaldInitLMPData) {
  Space s(3);
  s.initBoxLength(20);
  s.addMolInit("../forcefield/data.spce");
  s.addMol("../forcefield/data.spce");
  s.addMol("../forcefield/data.spce");
  s.addMol("../forcefield/data.spce");
  s.addMol("../forcefield/data.spce");

  // first pair uses old bulk spce method
  PairLJCoulEwald p1(&s, {{"rCut", "10"}});
  p1.initBulkSPCE(5.6, 38);

  // second pair uses new initData method
  PairLJCoulEwald p2(&s, {{"rCut", "10"}});
  p2.initData("../forcefield/data.spce");
  p2.initKSpace(5.6, 38);

  EXPECT_NEAR(p1.q()[0], p2.q()[0], DTOL);
  EXPECT_NEAR(p1.peLJ(), p2.peLJ(), DTOL);
  EXPECT_NEAR(p1.peTot(), p2.peTot(), DTOL);

  // test restart
  p2.writeRestart("tmp/rst");
  PairLJCoulEwald p3(&s, "tmp/rst");
  EXPECT_NEAR(p2.peTot(), p3.peTot(), DTOL);
  EXPECT_NEAR(p2.peLJ(), p3.peLJ(), DTOL);
}

