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
#include "mc_wltmmc.h"
#include "pair_hard_sphere.h"
#include "analyze_scatter.h"
#include "ui_abbreviated.h"
#include "trial_transform.h"

using namespace feasst;

void testVec(vector<vector<vector<vector<long long> > > > vec, vector<vector<vector<vector<long long> > > > vec2) {
  EXPECT_EQ(vec.size(), vec2.size());
  for (int i = 0; i < int(vec.size()); ++i) {
    EXPECT_EQ(vec[i].size(), vec2[i].size());
    for (int j = 0; j < int(vec[i].size()); ++j) {
      EXPECT_EQ(vec[i][j].size(), vec2[i][j].size());
      for (int k = 0; k < int(vec[i][j].size()); ++k) {
        EXPECT_EQ(vec[i][j][k].size(), vec2[i][j][k].size());
        for (int l = 0; l < int(vec[i][j][k].size()); ++l) {
          EXPECT_EQ(vec[i][j][k][l], vec2[i][j][k][l]);
        }
      }
    }
  }
}

TEST(Analyze, constructANDproduction) {
  Space s(3);
  s.initBoxLength(90);
  PairHardSphere p(&s);
  p.initData("../forcefield/data.cg4_mab");
  s.updateCells(p.rCutMaxAll());
  p.Forces();
  CriteriaMetropolis c(1, exp(-1));
  MC mc(&s, &p, &c);
  mc.weight = 1;
  transformTrial(&mc, "translate", 5);
  transformTrial(&mc, "rotate", 5);
  //initConfigBias(&mc, "../forcefield/data.cg4_mab");
  mc.nMolSeek(20, "../forcefield/data.cg4_mab", 1e9);
  mc.initRestart("tmp/anrst", 1e3);
  shared_ptr<AnalyzeScatter> scat = make_shared<AnalyzeScatter>(&p);
  scat->initSANS(0.5);
  scat->initFreq(1e2);
  scat->initFileName("tmp/iq");
  scat->initPrintFreq(1e3);
  scat->initProduction(0);
  mc.initAnalyze(scat);
  mc.runNumTrials(4*1e3);
  EXPECT_EQ(scat->production(), 0);
  scat->write();
  scat->writeRestart("tmp/hrst");
  mc.initProduction();
  EXPECT_EQ(scat->production(), 1);
  AnalyzeScatter scat2(&p, "tmp/hrst");
  EXPECT_EQ(scat2.production(), 0);
  scat2.writeRestart("tmp/hrst2");
  testVec(scat->histInter(), scat2.histInter());
  testVec(scat->histIntra(), scat2.histIntra());
}
 
TEST(AnalyzeScatter, args) {
  Space space;
  PairHardSphere pair(&space);
  try {
    AnalyzeScatter analyze(&pair, {{"/not/a/proper/arg", "error"}});
    CATCH_PHRASE("is not recognized");
  }
}
