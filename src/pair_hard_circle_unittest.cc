#include <gtest/gtest.h>
#include "pair_hard_circle.h"

TEST(PairHardCircle, hardCircle) {
  Space s(2, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(6,dim);
  s.addMolInit("../forcefield/data.atom");

  // add two particles with a certain distance separation
  vector<double> x(s.dimen(), 0.);
  s.xAdd = x;
  s.addMol("../forcefield/data.atom");
  x[0] = 0.99;
  s.xAdd = x;
  s.addMol("../forcefield/data.atom");

  const double dCircle = 1., rDep = 0.05;
  PairHardCircle p(&s, dCircle+2*rDep);
  p.initLMPData("../forcefield/data.atom");
  p.initRDep(rDep);
  p.rCutijset(0,0,dCircle+2*rDep);
  p.initEnergy();
  EXPECT_LT(1e200, p.peTot());

  x[0] = 0.01;
  s.transMol(1, x);
  p.initEnergy();
  EXPECT_NEAR(-3.926608784061230000, p.peTot(), 100*doubleTolerance);

  s.transMol(1, x);
  p.initEnergy();
  EXPECT_NEAR(-3.357306156792230000, p.peTot(), 100*doubleTolerance);

  x[0] = 0.085;
  s.transMol(1, x);
  p.initEnergy();
  EXPECT_NEAR(-0.044482468027066000, p.peTot(), 100*doubleTolerance);

  x[0] = 0.005;
  s.transMol(1, x);
  p.initEnergy();
  EXPECT_NEAR(0, p.peTot(), 1e-6);

  x[0] = 0.0000000001;
  s.transMol(1, x);
  p.initEnergy();
  EXPECT_NEAR(0, p.peTot(), 100*doubleTolerance);

  // test 3-body
  x[0] = 0.9;
  s.transMol(1, x);
  x[0] = 1;
  s.xAdd = x;
  s.addMol("../forcefield/data.atom");
  p.addPart();
  p.initEnergy();
  EXPECT_NEAR(2*-3.926608784061230000, p.peTot(), sqrt(doubleTolerance));

  EXPECT_EQ(1, p.checkEnergy(1e-8, 0));
  EXPECT_EQ(1, p.checkEnergy(1e-8, 1));

}

