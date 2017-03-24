#include <gtest/gtest.h>
#include "pair_tabular_1d.h"
#include "pair_lj_multi.h"

TEST(PairTabular1D, printReadTable) {
  const double rCut = 1.08;
  Space ss(3, 0);
  string addMolTypeA("../forcefield/data.lj");
  string addMolTypeB("../forcefield/data.ljb");
  PairLJMulti pp(&ss, rCut);
  ss.addMolInit(addMolTypeA.c_str());
  pp.initLMPData(addMolTypeA.c_str());
  ss.addMolInit(addMolTypeB.c_str());
  pp.initLMPData(addMolTypeB.c_str());
  pp.initExpType(4);
  pp.setLambdaij(0,0, 0);    //lambdaAA=0
  pp.setLambdaij(0,1, 0.5);  //lambdaAB=0.5
  pp.setLambdaij(1,1, 0);  //lambdaBB=0
  pp.lrcFlag = 0;
  pp.linearShift(1);    // cut and shift
  pp.initEnergy();
  pp.printTable("tmp/ptab",100,0.99);

  Space s(3, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(4,dim);
  PairTabular1D p(&s);
  s.addMolInit(addMolTypeA.c_str());
  p.initLMPData(addMolTypeA.c_str());
  s.addMolInit(addMolTypeB.c_str());
  p.initLMPData(addMolTypeB.c_str());
  p.readTable("tmp/ptab");
  p.initEnergy();
  s.writeRestart("tmp/rststab");
  p.writeRestart("tmp/rstptab");

  vector<double> x(s.dimen(), 0.);
  s.xAdd = x;
  s.addMol(addMolTypeA.c_str());
  x[1] = 1;
  s.xAdd = x;
  s.addMol(addMolTypeB.c_str());
  p.addPart();
}

