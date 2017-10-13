/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include <gtest/gtest.h>
#include "pair_hard_sphere.h"

using namespace feasst;

TEST(PairHardSphere, rCut) {
  Space space;
  space.lset(20.);
  space.addMolInit("../forcefield/data.lj");
  PairHardSphere pair(&space);
  pair.initData("../forcefield/data.lj");
  EXPECT_NEAR(pair.sig(0), 1., DTOL);
  EXPECT_NEAR(pair.eps(0), 1., DTOL);
  EXPECT_NEAR(pair.sigij(0, 0), 1., DTOL);
  EXPECT_NEAR(pair.epsij(0, 0), 1., DTOL);
  EXPECT_NEAR(pair.rCutij(0, 0), 1., DTOL);

  // put down two particles at known separation
  vector<double> xAdd(space.dimen());
  space.xAdd = xAdd;
  space.addMol();
  xAdd[0] = 1.05;
  space.xAdd = xAdd;
  space.addMol();
  pair.addPart();
  pair.initEnergy();
  EXPECT_NEAR(pair.peTot(), 0., DTOL);

  // move particle and test again
  space.xset(0.95, 1, 0);
  pair.initEnergy();
  EXPECT_NEAR(pair.peTot(), NUM_INF, DTOL);
}

