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
#include "pair_lj.h"
#include "pair_lj_coul_ewald.h"
#include "pair_patch_kf.h"
#include "pair_hybrid.h"

using namespace feasst;

TEST(PairHybrid, hybrid) {
  Space s;
  s.initBoxLength(24.8586887);
  s.readXYZBulk(3, "water", "../unittest/spce/test52.xyz");
  s.addMolInit("../forcefield/data.spce");
  PairLJCoulEwald p(&s, {{"rCut", "12.42934435"}});
  p.initBulkSPCE(5.6, 38);
  double petot = p.peTot();
  PairHybrid ph(&s);
  ph.addPair(&p);
  EXPECT_EQ(1, ph.nPairs());
  EXPECT_EQ(petot, ph.peTot());

  Space s2(3);
  s2.initBoxLength(8);
  s2.readXYZBulk(1, "atom", "../unittest/lj/srsw/lj_sample_config_periodic4.xyz");
  PairLJ p2(&s2, {{"rCut", "3"}, {"cutType", "linearShift"}});
  petot += p2.peTot();
  ph.addPair(&p2);
  EXPECT_EQ(2, ph.nPairs());
  EXPECT_EQ(petot, ph.peTot());

  Space s3(3);
  s3.initBoxLength(10);
  s3.readXYZBulk(2, "onePatch", "../unittest/patch/onePatch5.xyz");
  PairPatchKF p3(&s3, {{"rCut", "3"}, {"patchAngle", "90"}});
  p3.initEnergy();
  petot += p3.peTot();
  ph.addPair(&p3);
  EXPECT_EQ(3, ph.nPairs());
  EXPECT_EQ(petot, ph.peTot());

}

TEST(PairHybrid, ljANDwca) {

  Space s;
  s.initBoxLength(60);
  PairLJ pLJ(&s, {{"rCut", "3"}, {"molType", "../forcefield/data.cg3_60_43_1"}});
  vector<double> xAdd(3, 0.);
  pLJ.addMol(xAdd);
  xAdd[0] = 3.;
  pLJ.addMol(xAdd);
  xAdd[0] = -3.;
  pLJ.addMol(xAdd);
  xAdd[1] = 3.;
  pLJ.addMol(xAdd);
  xAdd[1] = -3.;
  pLJ.addMol(xAdd);
  pLJ.linearShift(1);
  pLJ.epsijset(1, 2, 0.);
  pLJ.epsijset(2, 1, 0.);
  pLJ.epsijset(2, 2, 0.);
  pLJ.initEnergy();
  PairLJ pWCA(&s, {{"rCut", feasst::str(pow(2, 1./6.))},
    {"molType", "../forcefield/data.cg3_60_43_1"}});
  pWCA.cutShift(1);
  pWCA.epsijset(1, 1, 0.);
  pWCA.initEnergy();
  PairHybrid ph(&s);
  ph.addPair(&pLJ);
  ph.addPair(&pWCA);
  ph.initEnergy();
  EXPECT_EQ(2, ph.nPairs());
  //EXPECT_NEAR(0.752213404790411000, ph.peTot(), 1e-16);
  EXPECT_EQ(1, pLJ.checkEnergy(1e-10, 1));
  EXPECT_EQ(1, pWCA.checkEnergy(1e-10, 1));
  EXPECT_EQ(1, ph.checkEnergy(1e-10, 1));
  EXPECT_EQ(1, pLJ.checkEnergy(1e-10, 0));

  // cout << "about to move" << endl;
  vector<int> mpart(4);
  mpart[0] = 0; mpart[1] = 1; mpart[2] = 2; mpart[3] = 3;
  ph.multiPartEner(mpart, 0);
  ph.update(mpart, 0, "store");
  s.randDisp(mpart, 2);
  s.randRotate(mpart, 2);
  ph.multiPartEner(mpart, 1);
  ph.update(mpart, 0, "update");
  // cout << "moved" << endl;
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
  PairLJ pLJ3(&s, "tmp/pLJrst");
  pWCA.writeRestart("tmp/pWCArst");
  PairLJ pWCA3(&s, "tmp/pWCArst");
  ph.writeRestart("tmp/phrst");
  PairHybrid ph3(&s, "tmp/phrst");

}

TEST(PairHybrid, Intra) {
  feasst::Space space;
  space.initBoxLength(90.);

  space.addMolInit("../forcefield/data.cg7mab2fullangflex");
  vector<double> xAdd(space.dimen(), 0.);
  space.xAdd = xAdd;
  space.addMol();
  xAdd[0] = 11;
  xAdd[1] = 5;
  xAdd[2] = 3;
  space.xAdd = xAdd;
  space.addMol();

  feasst::PairLJ pairInter(&space,
    {{"rCut", "30"},
     {"molType", "../forcefield/data.cg7mab2fullangflex"},
     {"cutType", "cutShift"}});
  feasst::PairLJ pairIntra(&space,
    {{"rCut", "30"},
     {"molType", "../forcefield/data.cg7mab2fullangflex"},
     {"cutType", "cutShift"}});

  // initialize intramolecular interactions via map
  vector<vector<int> > map(7, vector<int>(7, 1));
  for (int i = 0; i < 7; ++i) map[i][i] = 0;
  map[0][1] = 0;
  map[0][3] = 0;
  map[0][5] = 0;
  map[1][2] = 0;
  map[3][4] = 0;
  map[5][6] = 0;
  pairIntra.initIntra(2, map);

  feasst::PairHybrid pairHybrid(&space);
  pairHybrid.addPair(&pairInter);
  pairHybrid.addPair(&pairIntra);

  pairHybrid.initEnergy();
  EXPECT_NEAR(pairHybrid.peTot(), -7.3557992487162247, 10*feasst::DTOL);
  EXPECT_NEAR(pairInter.multiPartEner(space.imol2mpart(0), 0),
    pairInter.peTot(), feasst::DTOL);
  EXPECT_NEAR(pairIntra.multiPartEner(space.imol2mpart(0), 0) +
              pairIntra.multiPartEner(space.imol2mpart(1), 0),
    pairIntra.peTot(), 10*feasst::DTOL);
}

