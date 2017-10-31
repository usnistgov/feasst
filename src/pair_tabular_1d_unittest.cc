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
#include "pair_tabular_1d.h"
#include "pair_lj_multi.h"

using namespace feasst;

TEST(PairTabular1D, printReadTable) {
  const double rCut = 1.08;
  Space ss(3, 0);
  ss.initBoxLength(4);
  string addMolTypeA("../forcefield/data.lj");
  string addMolTypeB("../forcefield/data.ljb");
  PairLJMulti pp(&ss, rCut);
  PairTabular1D p(&ss);
  ss.addMolInit(addMolTypeA.c_str());
  pp.initData(addMolTypeA.c_str());
  p.initData(addMolTypeA.c_str());
  ss.addMolInit(addMolTypeB.c_str());
  pp.initData(addMolTypeB.c_str());
  p.initData(addMolTypeB.c_str());
  pp.initExpType(4);
  pp.setLambdaij(0,0, 0);    //lambdaAA=0
  pp.setLambdaij(0,1, 0.5);  //lambdaAB=0.5
  pp.setLambdaij(1,1, 0);  //lambdaBB=0
  pp.lrcFlag = 0;
  pp.linearShift(1);    // cut and shift
  pp.initEnergy();
  pp.printTable("tmp/ptab",100,0.99);
  p.readTable("tmp/ptab");
  p.initEnergy();

  ss.writeRestart("tmp/rststab");
  p.writeRestart("tmp/rstptab");

  vector<double> x(ss.dimen(), 0.);
  ss.xAdd = x;
  ss.addMol(addMolTypeA.c_str());
  x[1] = 1.05;
  ss.xAdd = x;
  ss.addMol(addMolTypeB.c_str());
  pp.addPart();
  p.addPart();

  pp.initEnergy();
  p.initEnergy();
  EXPECT_NEAR(pp.peTot(), -0.0033921341547825351, DTOL);
  EXPECT_NEAR(p.peTot(), -0.0033970032627577791, DTOL);
}

#ifdef GSL_
TEST(PairTabular1D, interpolateForces) {
  const double rCut = 3.;
  Space s(3, 0);
  s.initBoxLength(6);
  string addMolType("../forcefield/data.lj");
  PairLJMulti pp(&s, rCut);
  s.addMolInit(addMolType.c_str());
  pp.initData(addMolType.c_str());
  pp.cutShift(1);
  pp.initEnergy();
  pp.printTable("tmp/ptab",1000,0.99);

  PairTabular1D p(&s);
  p.initData(addMolType.c_str());
  p.readTable("tmp/ptab");
  p.setInterpolator("gslspline");
  p.initEnergy();

  vector<double> x(s.dimen(), 0.);
  s.xAdd = x;
  s.addMol(addMolType.c_str());
  x[1] = 1.2;
  s.xAdd = x;
  s.addMol(addMolType.c_str());
  p.addPart();
  p.initEnergy();
  pp.addPart();
  pp.initEnergy();

//  double fij = 0;
//  const double pe = p.peTable()[0][0]->interpolate(pow(x[1], 2), &fij);
//  cout << "pe " << pe << " fij " << fij << endl;
//  cout << "petab " << p.peTot() << " fijtab " << p.f(0,1) << endl;
//  cout << "pelj " << pp.peTot() << " fij " << pp.f(0,1) << endl;
  EXPECT_NEAR(p.peTot(), -0.88548584583883716, 1e-7);
  EXPECT_NEAR(p.f(0,1), 2.2116918590583787, 1e-5);
  EXPECT_NEAR(p.peTot(), pp.peTot(), 1e-7);
  EXPECT_NEAR(p.f(0,1), pp.f(0,1), 1e-5);
}
#endif  // GSL_

