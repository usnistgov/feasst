#include <gtest/gtest.h>
#include "mc.h"
#include "pair_hs.h"
#include "analyze_scatter.h"
#include "ui_abbreviated.h"

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

TEST(Analyze, construct) {
  Space s(3, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(90,dim);
  s.addMolInit("../forcefield/data.cg4_mab");
  s.initCellAtomCut(1);
  PairHS p(&s, 3.);
  p.initLMPData("../forcefield/data.cg4_mab");
  p.sig2rCut();
  s.updateCells(p.rCutMaxAll());
  p.Forces();
  CriteriaMetropolis c(1, exp(-1));
  MC mc(&s, &p, &c);
  mc.weight = 1;
  transformTrial(&mc, "translate", 5);
  transformTrial(&mc, "rotate", 5);
  initConfigBias(&mc, "../forcefield/data.cg4_mab");
  mc.nMolSeek(50, "../forcefield/data.cg4_mab", 1e9);
  mc.initRestart("tmp/anrst", 1e3);
  AnalyzeScatter scat(&s, &p);
  scat.initSANS(0.1);
  scat.initFreq(1e2);
  scat.initFileName("tmp/iq");
  scat.initPrintFreq(1e3);
  mc.initAnalyze(&scat);
  mc.runNumTrials(4*1e3);
  scat.print();
  scat.writeRestart("tmp/hrst");
  AnalyzeScatter scat2(&s, &p, "tmp/hrst");
  scat2.writeRestart("tmp/hrst2");
  testVec(scat.histInter(), scat2.histInter());
  testVec(scat.histIntra(), scat2.histIntra());
}

