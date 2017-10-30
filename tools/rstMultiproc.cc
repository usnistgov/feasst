/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "mc_wltmmc.h"

//const double newfCollect = 5e-5;
//
//class AnalyzeMonkeyPatch : public Analyze {
// public:
//  AnalyzeMonkeyPatch(Space *space, Pair *pair) : Analyze(space, pair) {}
//
//  void modifyRestart(shared_ptr<WLTMMC> mc) {
//    mc->c()->collectInit(newfCollect);
//    if (mc->c()->lnf() < newfCollect) {
//      mc->c()->collectInit();
//    }
//  }
//};

int main() {
  // set input variables
  std::ostringstream rstFileName("tmp/rst");

  // read restart file
  feasst::WLTMMC mc(rstFileName.str().c_str());

//  // patch
//  mc.c()->collectInit(newfCollect);
//  AnalyzeMonkeyPatch patch(mc.space(), mc.pair());
//  mc.initAnalyze(&patch);

  // run sweeps
  mc.runNumSweepsRestart(100, rstFileName.str().c_str());
}
