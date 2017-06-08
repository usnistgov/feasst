#include <gtest/gtest.h>
#include "mc_wltmmc.h"
#include "pair_hs.h"
#include "analyze_scatter.h"
#include "ui_abbreviated.h"

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
  Space s(3, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(90,dim);
  s.addMolInit("../forcefield/data.cg4_mab");
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
  //initConfigBias(&mc, "../forcefield/data.cg4_mab");
  mc.nMolSeek(20, "../forcefield/data.cg4_mab", 1e9);
  mc.initRestart("tmp/anrst", 1e3);
  AnalyzeScatter scat(&s, &p);
  scat.initSANS(0.5);
  scat.initFreq(1e2);
  scat.initFileName("tmp/iq");
  scat.initPrintFreq(1e3);
  scat.initProduction(0);
  mc.initAnalyze(&scat);
  mc.runNumTrials(4*1e3);
  EXPECT_EQ(scat.production(), 0);
  scat.write();
  scat.writeRestart("tmp/hrst");
  mc.initProduction();
  EXPECT_EQ(scat.production(), 1);
  AnalyzeScatter scat2(&s, &p, "tmp/hrst");
  EXPECT_EQ(scat2.production(), 0);
  scat2.writeRestart("tmp/hrst2");
  testVec(scat.histInter(), scat2.histInter());
  testVec(scat.histIntra(), scat2.histIntra());
}

const double newfCollect = 5e-5;
class AnalyzeMonkeyPatch : public Analyze {
 public:
  //AnalyzeMonkeyPatch() : Analyze() {}
  AnalyzeMonkeyPatch(Space *space, Pair *pair) : Analyze(space, pair) {}
  void modifyRestart(shared_ptr<WLTMMC> mc) {
    mc->c()->collectInit(newfCollect);
    if (mc->c()->lnf() < newfCollect) {
      mc->c()->collectInit();
    }
  }
};

TEST(Analyze, MonkeyPatch) {
  Space s(3, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(90,dim);
  s.addMolInit("../forcefield/data.cg4_mab");
  PairHS p(&s, 3.);
  p.initLMPData("../forcefield/data.cg4_mab");
  p.sig2rCut();
  s.updateCells(p.rCutMaxAll());
  p.Forces();
  CriteriaWLTMMC c(1, exp(-1), "nmol", -0.5, 20.5, 20);
  WLTMMC mc(&s, &p, &c);
  mc.weight = 1;
  transformTrial(&mc, "translate", 5);
  transformTrial(&mc, "rotate", 5);
  mc.nMolSeek(20, "../forcefield/data.cg4_mab", 1e9);
  mc.initRestart("tmp/mpanrst", 1e1);
  mc.initColMat("tmp/monkeycol", 1e1);
  mc.initWindows(1);
  mc.writeRestart("tmp/monkeyrst");
  int t = 0;
  #ifdef _OPENMP
    #pragma omp parallel private(t)
    {
      t = omp_get_thread_num();
      stringstream ss;
      ss << "tmp/monkeyrstp" << t;
      mc.writeRestart(ss.str().c_str());
    }
    #pragma omp barrier

    WLTMMC mc2("tmp/monkeyrst"); 
    mc2.c()->collectInit(newfCollect);
    AnalyzeMonkeyPatch patch(mc2.space(), mc2.pair());

    mc2.runNumSweepsRestart(0, "tmp/monkeyrst");

    CriteriaWLTMMC c2("tmp/monkeyrstp0criteria");
    EXPECT_NEAR(1e-6, c2.lnfCollect(), DTOL);
    EXPECT_NE(newfCollect, c2.lnfCollect());
    
    // apply patch
    mc2.initAnalyze(&patch);
    mc2.runNumSweepsRestart(0, "tmp/monkeyrst");
    CriteriaWLTMMC c3("tmp/monkeyrstp0criteria");
    EXPECT_NEAR(newfCollect, c3.lnfCollect(), DTOL);
  #endif  // _OPENMP
}

