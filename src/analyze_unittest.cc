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
#include "analyze.h"
#include "trial_transform.h"
#include "ui_abbreviated.h"

using namespace feasst;

const double newfCollect = 5e-5;
class AnalyzeMonkeyPatch : public Analyze {
 public:
  //AnalyzeMonkeyPatch() : Analyze() {}
  AnalyzeMonkeyPatch(Pair *pair) : Analyze(pair) {}
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
  PairHardSphere p(&s);
  p.initData("../forcefield/data.cg4_mab");
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
  #ifdef _OPENMP
    int t = 0;
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
    shared_ptr<AnalyzeMonkeyPatch> patch = make_shared<AnalyzeMonkeyPatch>(mc2.pair());

    mc2.runNumSweepsRestart(0, "tmp/monkeyrst");

    CriteriaWLTMMC c2("tmp/monkeyrstp0criteria");
    EXPECT_NEAR(1e-6, c2.lnfCollect(), DTOL);
    EXPECT_NE(newfCollect, c2.lnfCollect());

    // apply patch
    mc2.initAnalyze(patch);
    mc2.runNumSweepsRestart(0, "tmp/monkeyrst");
    CriteriaWLTMMC c3("tmp/monkeyrstp0criteria");
    EXPECT_NEAR(newfCollect, c3.lnfCollect(), DTOL);
  #endif  // _OPENMP
}
