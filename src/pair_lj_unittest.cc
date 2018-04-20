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

using namespace feasst;

TEST(PairLJ, dim) {
  for (int dimen = 2; dimen < 4; ++dimen) {
    for (int natom = 1; natom != 5; ++natom) {
      Space s(dimen);
      s.init_config(natom);
      PairLJ p(&s);
      vector<double> pe = p.pe();
      EXPECT_EQ(natom,int(pe.size()));
      vector<vector<double> > f = p.f();
      EXPECT_EQ(natom,int(f.size()));
      vector<vector<vector<double> > > vr = p.vr();
      EXPECT_EQ(natom,int(vr.size()));
      for (int i = 0; i < natom; ++i) {
        EXPECT_EQ(dimen,int(f[i].size()));
        EXPECT_EQ(dimen,int(vr[i].size()));
        for (int j = 0; j < dimen; ++j) {
          EXPECT_EQ(dimen,int(vr[i][j].size()));
        }
      }
    }
  }
}

// Check that total potential energy using PairLJ gives same result as multiPartEner calculated for all particles individually, and halved

TEST(PairLJ, checkEnergy) {
  Space s(3);
  s.init_config(12);
  feasst::PairLJ p(&s, {{"rCut", "5"}, {"cutType", "none"}});
  EXPECT_EQ(1, p.checkEnergy(1e-11, 1));
}

TEST(PairLJ, addPartdelPart) {
  const int dim = 3, natom = 12;
  Space s(dim);
  s.init_config(natom);
  PairLJ p(&s, {{"rCut", "5"}, {"cutType", "none"}});

  // remove particle
  vector<int> mpart(1, 0);
  p.delPart(mpart);
  s.delPart(mpart);
  p.initEnergy();
  EXPECT_EQ(11, int(p.f().size()));
  EXPECT_EQ(11, int(p.vr().size()));
  EXPECT_EQ(11, int(p.pe().size()));

  // add particle
  vector<double> x(s.dimen());
  x.at(0) = 0; x.at(1) = 0; x.at(2) = 0;
  s.addPart(x,0,1);
  p.addPart();
  p.initEnergy();
  EXPECT_EQ(12, int(p.f().size()));
  EXPECT_EQ(12, int(p.vr().size()));
  EXPECT_EQ(12, int(p.pe().size()));

}

TEST(PairLJ, equltl43muvttmmc) {
  Space s(3);
  s.initBoxLength(9);
  s.readXYZBulk(4, "equltl43", "../unittest/equltl43/two.xyz");
  PairLJ p(&s, {{"rCut", feasst::str(pow(2, 1./6.))}, {"molType", "../forcefield/data.equltl43"}});
  // energy should be equal to 3 wca beads touching at distance 1.1sig
  p.initWCA(1, 1);
  p.initEnergy();
  EXPECT_NEAR(3*(4*(pow(1.1,-12)-pow(1.1,-6))+1), p.peTot(), 1e-15);
  EXPECT_EQ(1, p.checkEnergy(1e-18, 1));
}

TEST(PairLJ, linearForceShiftLJ) {
  Space s(3);
  std::stringstream addMolType;
  addMolType << s.install_dir() << "/forcefield/data.lj";
  s.addMolInit(addMolType.str());
  std::ifstream file("../unittest/lj/srsw/lj_sample_config_periodic4.xyz");
  s.readXYZ(file);
  EXPECT_EQ(30, s.nMol());
  PairLJ p(&s, {
    {"rCut", "3"},
    {"cutType", "linearShift"},
    {"molType", addMolType.str()}});
  EXPECT_EQ(3, p.rCutij(0, 0));
  if (p.peTot() < 100) EXPECT_EQ(1, p.checkEnergy(1e-10, 1));
}

TEST(PairLJ, exVol) {
  Space s(3);
  s.initBoxLength(0);
  PairLJ p(&s);
  p.addMol();
  const double boxl = 2.*(s.maxMolDist() + 1 + 0.1);
  s.initBoxLength(boxl);
  EXPECT_EQ(0, s.x(0,0));
  EXPECT_NEAR(4.*PI/3., p.exVol(1e2), 3e-3);
}

TEST(PairLJ, args) {
  {
    feasst::Space space;
    feasst::PairLJ pair(&space, {{"rCut", "3."}});
    std::stringstream ss;
    ss << space.install_dir() << "/forcefield/data.lj";
    EXPECT_EQ(space.addMolListType().begin()->compare(ss.str()), 0);
    EXPECT_EQ(pair.lrcFlag, 1);
  }

  {
    feasst::Space space;
    std::stringstream ss;
    ss << space.install_dir() << "/forcefield/data.ljb";
    feasst::PairLJ pair(&space, {{"rCut", "3."}, {"cutType", "cutShift"},
      {"molType", ss.str()}});
    EXPECT_EQ(space.addMolListType()[0], ss.str());
    EXPECT_EQ(pair.lrcFlag, 0);
  }
  try {
    feasst::Space space;
    feasst::PairLJ p(&space, {{"not/a/proper/arg", "error"}});
    CATCH_PHRASE("is not recognized");
  }

  // cannot use "none" for molType with other arguments as well
  try {
    feasst::Space space;
    feasst::PairLJ p(&space, {{"molType", "none"}, {"cutType", "lrc"}});
    CATCH_PHRASE("is not recognized");
  }
}

TEST(PairLJ, WCAanalytical) {

  // WCA for sig=1 and 0.85
  for (double sig = 0.85; sig < 1.01; sig += 0.15) {
    Space s(3);
    s.initBoxLength(100);
    stringstream ss;
    ss << "../forcefield/data.lj";
    if (sig == 0.85) ss << "s0.85";
    PairLJ p(&s, {{"molType", ss.str()}});
    p.initWCA(0,0);
    vector<double> x(s.dimen(), 0.);
    p.addMol(x);
    x[0] = 1.;
    p.addMol(x);
    EXPECT_NEAR(1, s.x(1,0) - s.x(0,0), 1e-15);
    EXPECT_NEAR(0, s.x(1,1) - s.x(0,1), 1e-15);
    EXPECT_NEAR(0, s.x(1,2) - s.x(0,2), 1e-15);


    if (sig == 1) {
      s.xset(1, 1, 0);
      p.initEnergy();
      EXPECT_NEAR(1, p.peTot(), DTOL);
      s.xset(pow(2, 1./6.), 1, 0);
      p.initEnergy();
      EXPECT_NEAR(0, p.peTot(), DTOL);
      s.xset(1.01, 1, 0);
      p.initEnergy();
      EXPECT_NEAR(0.781615960043788000, p.peTot(), DTOL);
    } else {
      s.xset(0.85, 1, 0);
      p.initEnergy();
      EXPECT_NEAR(1, p.peTot(), DTOL);
      s.xset(pow(2, 1./6.)*0.85, 1, 0);
      p.initEnergy();
      EXPECT_NEAR(0, p.peTot(), DTOL);
      s.xset(0.9, 1, 0);
      p.initEnergy();
      EXPECT_NEAR(0.175851657413201000, p.peTot(), DTOL);
    }
  }
}

TEST(PairLJ, cg3analytical) {
  Space s(3);
  PairLJ p(&s, {{"rCut", "3"}, {"molType", "../forcefield/data.cg3_60_1_1"}});
  vector<double> x(s.dimen(), 0.);
  p.addMol(x);
  x[0] = 2.;
  p.addMol(x);
  p.rCutijset(1, 1, p.rCut());
  p.linearShiftijset(1, 1, 1);
  p.initWCA(1, 2);
  p.initWCA(2, 2);
  p.initEnergy();
  EXPECT_NEAR(p.peTot(), 0.95489983432133096, DTOL);
}

TEST(PairLJ, LJYanalytical) {
  Space s(3);
  PairLJ p(&s, {{"rCut", "3"}});
  vector<double> x(s.dimen(), 0.);
  p.addMol(x);
  x[0] = 2.;
  p.addMol(x);
  p.initExpType(1);
  p.initScreenedElectro(2, 0.5);
  p.linearShift(0);
  p.lrcFlag = 0;
  p.initEnergy();
  EXPECT_NEAR(p.peTot(), 0.366903117090021000, DTOL);
  p.linearShift(1);
  p.initEnergy();
  EXPECT_NEAR(p.peTot(), 0.094226110309604982, DTOL);
}

TEST(PairLJ, MMLJanalytical) {
  Space s(3);
  PairLJ p(&s, {{"rCut", "3"}, {"molType", "../forcefield/data.cg3_91_0.57_2"}});
  vector<double> x(s.dimen(), 0.);
  p.addMol(x);
  x[0] = 2.*0.85;
  p.addMol(x);
  p.rCutijset(1, 1, 3);
  p.linearShiftijset(1, 1, 1);
  p.initWCA(1, 2);
  p.initWCA(2, 2);
  p.initEnergy();
  EXPECT_NEAR(p.peTot(), -0.139144838318008000 + 0.302421948569261000, 5*DTOL);
  x[0] = -0.9;
  x[1] = 1;
  s.transMol(1, x);
  p.initEnergy();
  EXPECT_NEAR(p.peTot(), -0.676951777428172000 + 66.449659060736600000, 1000*DTOL);
}

TEST(PairLJ, cg3exampleConfig) {
  Space s(3);
  s.addMolInit("../forcefield/data.cg3_60_1_1");
  std::ifstream inFile("../unittest/cg3/cg3_60_1_1/example/moviep1n50.xyz");
  s.readXYZ(inFile);
  EXPECT_EQ(50, s.nMol());
  vector<double> x(s.dimen(), 0.);
  PairLJ p(&s, {{"rCut", "3"}, {"molType", "../forcefield/data.cg3_60_1_1"}});
  p.rCutijset(1, 1, p.rCut());
  p.linearShiftijset(1, 1, 1);
  p.rCutijset(1, 2, pow(2, 1./6.));
  p.cutShiftijset(1, 2, 1);
  p.rCutijset(2, 2, pow(2, 1./6.));
  p.cutShiftijset(2, 2, 1);
  p.initEnergy();
  EXPECT_NEAR(p.peTot(), -108.25895014589899, DTOL);
}

TEST(PairLJ, cg3analyticalAlpha128) {
  Space s(3);
  PairLJ p(&s, {{"rCut", "1.08"}, {"molType", "../forcefield/data.cg3_91_0.57_2"}});
  vector<double> x(s.dimen(), 0.);
  p.addMol(x);
  x[1] = 1.22;
  p.addMol(x);

  p.initExpType(4); // alpha=128
  p.rCutijset(1, 1, p.rCut());
  p.linearShiftijset(1, 1, 1);
  p.initWCA(1, 2);
  p.initWCA(2, 2);
  p.initEnergy();

  EXPECT_NEAR(p.peTot(), 2*39.789289254425900000, 10000*DTOL);
  //EXPECT_NEAR(p.f(0, 1), 0, DTOL);
  //EXPECT_NEAR(sqrt(pow(p.f(1, 1),2)+pow(p.f(1, 0), 2)), 23095.254845213545, 1e-9);

  // test peMap and neighCut
  const double peTot = p.peTot();
  p.initNeighCut(1);
  p.initPEMap(1);
  const vector<int> mpart = s.randMol();
  p.multiPartEner(mpart, 0);
  vector<int> neigh;
  vector<double> peMap;
  p.neighCutMolPEMap(neigh, peMap);
  EXPECT_EQ(1, int(neigh.size()));
  EXPECT_EQ(1, int(peMap.size()));
  EXPECT_NEAR(peMap[0], peTot, DTOL);

  // flip
  s.qMolAlt(1, 0, 1);
  s.qMolAlt(1, 3, 0);
  s.quat2pos(1);
  x[1] = -s.x(4, 1) + 2*0.266345520433943000 + 1.02;
  s.transMol(1, x);
  p.initEnergy();
  EXPECT_NEAR(p.peTot(), -0.290278106387070000, 100*DTOL);
  //EXPECT_NEAR(p.f(0, 1), 0, DTOL);
  //EXPECT_NEAR(sqrt(pow(p.f(1, 1),2)+pow(p.f(1, 0), 2)), 33.461405536957100000, 20000*DTOL);
  EXPECT_EQ(1, p.checkEnergy(DTOL, 2));
}

TEST(PairLJ, Gaussian) {
  Space s(3);
  PairLJ p(&s, {{"rCut", "3."}});
  vector<double> x(s.dimen(), 0.);
  p.addMol(x);
  x[1] = pow(2,1./6.);
  p.addMol(x);

  p.lrcFlag = 0;
  p.cutShift(0);
  p.initEnergy();
  EXPECT_NEAR(-1, p.peTot(), 10000*DTOL);

  const double height = 1.3408572374689746728347876674875;
  p.addGaussian(height, x[1], 1);

  p.initEnergy();
  EXPECT_NEAR(-1+height, p.peTot(), 10000*DTOL);

  p.addGaussian(1, 1, 2);
  p.initEnergy();
  EXPECT_NEAR(-1+height+exp(-pow(((x[1]-1)/2),2)), p.peTot(), 10000*DTOL);

  p.writeRestart("tmp/yoyorst");
  PairLJ p2(&s, "tmp/yoyorst");
  p2.writeRestart("tmp/yoyorst2");
  p2.initEnergy();

  EXPECT_NEAR(p.peTot(), p2.peTot(), 10000*DTOL);
}

TEST(PairLJ, Lambda) {
  Space s(2);
  s.initBoxLength(9);
  PairLJ p(&s, {{"rCut", "3."}});
  vector<double> x(s.dimen(), 0.);
  p.addMol(x);
  x[1] = pow(2,1./6.);
  p.addMol(x);

  p.setLambdaij(0,0,1);
  p.lrcFlag = 0;
  p.cutShift(1);
//  p.cutShift(0);
//  p.rCutijset(0,0,rCut);
//  p.initLRC();
  p.initEnergy();
  EXPECT_NEAR(-1+0.005479441744238780, p.peTot(), 10000*DTOL);
  p.setLambdaij(0,0,-1);
  p.initEnergy();
  EXPECT_NEAR(1-0.005479441744238780, p.peTot(), 10000*DTOL);
  p.setLambdaij(0, 0, 0);
  p.initEnergy();
  EXPECT_NEAR(0., p.peTot(), 10000*DTOL);

//  exit(0);
  p.writeRestart("tmp/yoyorst");
  PairLJ p2(&s, "tmp/yoyorst");
  p2.writeRestart("tmp/yoyorst2");
  p2.initEnergy();

  EXPECT_NEAR(p.peTot(), p2.peTot(), 10000*DTOL);
}

TEST(PairLJ, sigrefAnalytical) {
  Space s(3);
  PairLJ p(&s, {{"rCut", "3.3"}, {"molType", "../forcefield/data.ljs0.85"}});
  p.setSigRefFlag(1);
  p.initData("../forcefield/data.ljs0.85");

  vector<double> x(s.dimen(), 0.);
  p.addMol(x);
  x[1] = 1.22;
  p.addMol(x);

  p.initExpType(1); // alpha1:12
  p.setLambdaij(0, 0, -1);

  PairLJ *pcut = p.clone(&s);

  p.linearShift(1);
  //p.cutShift(1);
  p.initEnergy();

  EXPECT_NEAR(p.sigRef(0), 2.1, DTOL);
  EXPECT_NEAR(p.peTot(), 0.4867769858755370, 10000*DTOL);

  p.setLambdaij(0, 0, 1);
  p.initEnergy();
  EXPECT_NEAR(p.peTot(), -0.4867769858755370, 10000*DTOL);

  p.setLambdaij(0, 0, 0.331);
  p.initEnergy();
  EXPECT_NEAR(p.peTot(), -0.1611231823248030, 10000*DTOL);

  x[1] = -s.x(1, 1) + 2*x[1];
  s.transMol(1, x);
  p.initEnergy();
  EXPECT_NEAR(p.peTot(), -0.0011223005336472365, 10000*DTOL);

  x[1] = -s.x(1, 1) + 0.87;
  s.transMol(1, x);
  // test linear shift
  p.setLambdaij(0,0,-1);
  p.initEnergy();
  EXPECT_NEAR(p.peTot(), 1.6134113229466700, 10000*DTOL);
  // test cut shift
  pcut->setLambdaij(0,0,-1);
  pcut->cutShift(1);
  pcut->initEnergy();
  EXPECT_NEAR(pcut->peTot(), 1.6158060153131700, 10000*DTOL);

  x[1] = -s.x(1, 1) + 0.5;
  s.transMol(1, x);
  p.initEnergy();
  EXPECT_NEAR(p.peTot(), 284.3198540118970000, 10000*DTOL);

  delete pcut;
}

TEST(PairLJ, InitAlphaAnalytical) {
  Space s(3);
  PairLJ p(&s, {{"rCut", "1e4"}});

  vector<double> x(s.dimen(), 0.);
  p.addMol(x);
  x[1] = 1.22;
  p.addMol(x);

  p.initAlpha(5.5);

  PairLJ *pcut = p.clone(&s);
  p.linearShift(1);
  //p.cutShift(1);
  p.initEnergy();

  EXPECT_NEAR(p.peTot(), -0.8910756889503104, 10*DTOL);

  x[1] = -s.x(1, 1) + 2*x[1];
  s.transMol(1, x);
  p.initEnergy();
  EXPECT_NEAR(p.peTot(), -0.029389303309848992, 10*DTOL);

  delete pcut;
}
