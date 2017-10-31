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
#include "pair_hard_circle.h"

using namespace feasst;

TEST(PairHardCircle, hardCircle) {
  Space s(2);
  s.initBoxLength(6);
  const double dCircle = 1., rDep = 0.05;
  PairHardCircle p(&s, dCircle+2*rDep);
  p.initData("../forcefield/data.atom");

  // add two particles with a certain distance separation
  vector<double> x(s.dimen(), 0.);
  p.addMol(x);
  x[0] = 0.99;
  p.addMol(x);

  p.initRDep(rDep);
  p.rCutijset(0,0,dCircle+2*rDep);
  p.initEnergy();
  EXPECT_LT(1e200, p.peTot());

  x[0] = 0.01;
  s.transMol(1, x);
  p.initEnergy();
  EXPECT_NEAR(-3.926608784061230000, p.peTot(), 100*DTOL);

  s.transMol(1, x);
  p.initEnergy();
  EXPECT_NEAR(-3.357306156792230000, p.peTot(), 100*DTOL);

  x[0] = 0.085;
  s.transMol(1, x);
  p.initEnergy();
  EXPECT_NEAR(-0.044482468027066000, p.peTot(), 100*DTOL);

  x[0] = 0.005;
  s.transMol(1, x);
  p.initEnergy();
  EXPECT_NEAR(0, p.peTot(), 1e-6);

  x[0] = 0.0000000001;
  s.transMol(1, x);
  p.initEnergy();
  EXPECT_NEAR(0, p.peTot(), 100*DTOL);

  // test 3-body
  x[0] = 0.9;
  s.transMol(1, x);
  x[0] = 1;
  p.addMol(x);
  p.initEnergy();
  EXPECT_NEAR(2*-3.926608784061230000, p.peTot(), sqrt(DTOL));

  EXPECT_EQ(1, p.checkEnergy(1e-8, 0));
  EXPECT_EQ(1, p.checkEnergy(1e-8, 1));

}

