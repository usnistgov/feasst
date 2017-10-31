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
#include <limits.h>
#include "criteria_mayer.h"
#include "pair_lj_multi.h"
#include "pair_hard_sphere.h"
#include "pair_squarewell.h"
#include "trial_transform.h"
#include "mc.h"

TEST(CriteriaMayer, ref) {
  feasst::Space space(3);
  space.initBoxLength(8);
  feasst::PairSquareWell pair(&space, 2.);
  pair.initData("../forcefield/data.lj");
  vector<double> xAdd(3, 0.);
  pair.addMol(xAdd);
  xAdd[0] = 1.1;
  pair.addMol(xAdd);
  pair.rCutijset(0, 0, 2.);
  pair.initEnergy();

  const double beta = 1.;
  feasst::CriteriaMayer crit(beta);
  crit.initPairRef(&pair);
  crit.store(&pair);
  EXPECT_EQ(0, crit.accept(0, 0., "move", 0));

}

//TEST(MC, b2mayer) {
//  feasst::ranInitByDate();
//  feasst::Space space(3, 0);
//  space.addMolInit("../forcefield/data.lj");
//  vector<double> xAdd(space.dimen(), 0.);
//  space.xAdd = xAdd;
//  space.addMol(0);
//  xAdd[0] = 1.1;
//  space.xAdd = xAdd;
//  space.addMol(0);
//
//  feasst::PairLJMulti pair(&space, 1e5);
//  pair.initData("../forcefield/data.lj");
//  pair.cutShift(1);
//  pair.initEnergy();
//  feasst::PairHardSphere pairRef(&space, 1);
//  pairRef.initData("../forcefield/data.lj");
//  pairRef.initEnergy();
//
//  const double boxl = 2.*(2.*space.maxMolDist() + pair.rCut());
//  space.initBoxLength(boxl);
//
//  feasst::CriteriaMayer crit(0.2);
//  crit.initPairRef(&pairRef);
//  feasst::MC mc(&space, &pair, &crit);
//  transformTrial(&mc, "translate", 0.1);
//  //mc.setNFreqTune(1e5);
//  const int nfreq = 1e0;
//  mc.setNFreqCheckE(1, 1e10);
//  //mc.setNFreqCheckE(nfreq, 1e-8);
//  mc.initLog("tmp/mayer", nfreq);
//  mc.initMovie("tmp/mayer", nfreq);
//
//  //mc.nMolSeek(2, "../forcefield/data.lj", 1e8);
//
//  mc.runNumTrials(1e4);
//  //mc.runNumTrials(1e6);
//  //cout << "pe " << pair.peTot() << endl;
//  EXPECT_NEAR(crit.b2ratio(), 0.2, 1.);
//  // beta = 1 takes too long
//  // // EXPECT_NEAR(crit.b2ratio(), -2.5381, DTOL);
//}
//
