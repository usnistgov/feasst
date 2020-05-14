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
#include "pair_squarewell.h"

using namespace feasst;

TEST(PairSquareWell, mickeymouse) {
  Space s(3);
  s.initBoxLength(20);
  PairSquareWell p(&s, {{"rCut", "1.02"}});
  p.initData("../forcefield/data.cg3_91_0.57_2");
  vector<double> xAdd(s.dimen());
  p.addMol(xAdd);
  xAdd[1] = 2*0.266345520433943000 + 1.01;
  p.addMol(xAdd);
  // flip
  s.qMolAlt(1, 0, 1);
  s.qMolAlt(1, 3, 0);
  s.quat2pos(1);
  p.rCutijset(1, 1, p.rCut());
  p.initHardSphere(1, 2);
  p.initHardSphere(2, 2);
  p.initEnergy();
  EXPECT_NEAR(p.peTot(), -1, DTOL);
  EXPECT_EQ(p.checkEnergy(DTOL, 0), 1);

  // write restart, read restart and test
  { p.writeRestart("tmp/p");
    PairSquareWell p2(&s, "tmp/p");
    p2.initEnergy();
    EXPECT_NEAR(p2.peTot(), -1, DTOL);
  }

  xAdd[1] = 0.01000001;
  s.transMol(1, xAdd);
  p.initEnergy();
  EXPECT_NEAR(p.peTot(), 0, DTOL);
  EXPECT_EQ(p.checkEnergy(DTOL, 1), 1);
  xAdd[1] = -2*0.01000001;
  s.transMol(1, xAdd);
  p.initEnergy();
  EXPECT_GE(p.peTot(), 1e100);

  // flip particle again
  s.qMolAlt(1, 0, 0);
  s.qMolAlt(1, 3, 1);
  s.quat2pos(1);
  xAdd[1] = -s.x(4, 1) + 0.83 + 1*0.3966345520433943000;
  s.transMol(1, xAdd);
  p.initEnergy();
  EXPECT_GE(p.peTot(), 1e100);
  xAdd[1] = 0.02;
  s.transMol(1, xAdd);
  p.initEnergy();
  EXPECT_EQ(p.peTot(), 0);
  EXPECT_EQ(p.checkEnergy(DTOL, 1), 1);

  // write restart, read restart and test
  { p.writeRestart("tmp/p");
    PairSquareWell p2(&s, "tmp/p");
    p2.initEnergy();
    EXPECT_EQ(p2.peTot(), 0);
  }

  // test again
  xAdd[0] = s.x(2,0);
  xAdd[1] = -s.x(4,1)-0.133172760216972000-0.5*1.85-s.x(1,1)+1e-10;
  s.transMol(1, xAdd);
  p.initEnergy();
  EXPECT_GE(p.peTot(), 1e100);
  xAdd[0] = 0.;
  xAdd[1] = -2e-10;
  s.transMol(1, xAdd);
  p.initEnergy();
  EXPECT_EQ(p.peTot(), 0);
}

