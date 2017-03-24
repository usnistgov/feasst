#include <gtest/gtest.h>
#include "pair_lj.h"

TEST(PairLJ, dim) {
  for (int dimen = 1; dimen != 4; ++dimen) {
    for (int natom = 1; natom != 5; ++natom) {
      Space s(dimen,0);
      s.init_config(natom);
      PairLJ p(&s, 0.01);
      p.initEnergy();
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

// Check Forces between two particles against analytical expression

TEST(PairLJ, analytical) {
  const double boxl = 100;  // large domain size
  const double rCut = 4.;
  for (int dimen = 1; dimen != 2; ++dimen) {
    Space s(dimen,0);
    s.init_config(2);
    for (int dim = 0; dim < s.dimen(); ++dim) s.lset(boxl, dim);
    PairLJ p(&s, rCut);
    p.lrcFlag = 0;
    for (int i = 0; i < dimen; ++i) {
      s.lset(boxl, i);
    }
    s.xset(49., 0, 0);
    s.xset(47., 1, 0);
    while (s.x()[dimen*1+0] > -50.) {
      p.initEnergy();
      if (s.x()[dimen*1+0] == -1.) {
        // particle outside of cuttoff
        EXPECT_EQ(0, p.pe()[1]);
        EXPECT_EQ(0, p.f()[1][0]);
      } else if (s.x(1,0) == 49.) {
        // particle right next to other
//        EXPECT_NEAR(-0.0605471, pair->pe()[1]*2., 1e-7);
        EXPECT_NEAR(-0.061523455, p.pe()[1]*2., 1e-7);
        EXPECT_NEAR(0.1816406, p.f()[1][0], 1e-7);
      } else if (s.x(1,0) == -49.) {
        // particle right next to other, across pbc
        EXPECT_NEAR(-0.061523455, p.pe()[1]*2., 1e-7);
        EXPECT_NEAR(-0.1816406, p.f()[1][0], 1e-7);
      }
      s.xset(s.x(1,0) - 48., 1, 0);
    }
  }
}

// Check that total potential energy using PairLJ gives same result as multiPartEner calculated for all particles individually, and halved

TEST(PairLJ, checkEnergy) {
  Space s(3,0);
  s.init_config(12);
  PairLJ p(&s, 5);
  p.lrcFlag = 0;
  EXPECT_EQ(1, p.checkEnergy(1e-11, 1));
}

TEST(PairLJ, addPartdelPart) {
  const int dim=3, natom=12, rCut=5;
  Space s(dim,0);
  s.init_config(natom);
  PairLJ p(&s, rCut);

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
  Space s(3, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(9,dim);
  s.readXYZBulk(4, "equltl43", "../unittest/equltl43/two.xyz");
  PairLJ p(&s, pow(2, 1./6.));
  p.cutShift(1);
  p.lrcFlag = 0;
  p.initLMPData("../forcefield/data.equltl43");
  p.initEnergy();
  // energy should be equal to 3 wca beads touching at distance 1.1sig
  EXPECT_NEAR(3*(4*(pow(1.1,-12)-pow(1.1,-6))+1), p.peTot(), 1e-15);
  EXPECT_EQ(1, p.checkEnergy(1e-18, 1));
}

TEST(PairLJ, linearForceShiftLJ) {
  Space s(3, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(8,dim);
  s.readXYZBulk(1, "atom", "../unittest/lj/srsw/lj_sample_config_periodic4.xyz");
  s.addMolInit("../forcefield/data.lj");
  PairLJ p(&s, 3);
  p.linearShift(1);
  p.initEnergy();
  s.addMol("../forcefield/data.lj");
  p.addPart();
  p.initEnergy();
  if (p.peTot() < 100) EXPECT_EQ(1, p.checkEnergy(1e-10, 1));
}

TEST(PairLJ, exVol) {
  Space s(3, 0);
  for (int dim = 0; dim < s.dimen(); ++dim) s.lset(0, dim);
  s.addMolInit("../forcefield/data.lj");
  s.addMol("../forcefield/data.lj");
  const double boxl = 2.*(s.maxMolDist() + 1 + 0.1);
  for (int dim = 0; dim < s.dimen(); ++dim) s.lset(boxl, dim);
  EXPECT_EQ(0, s.x(0,0));
  PairLJ p(&s, 1e-12);
  p.initLMPData("../forcefield/data.lj");
  EXPECT_NEAR(4.*PI/3., p.exVol(1e2), 3e-3);
  //EXPECT_NEAR(4.*PI/3., p.exVol(1e3), 2e-4);
}


