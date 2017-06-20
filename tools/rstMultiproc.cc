/**
 * \mainpage
 *
 * Developed by Harold Wickes Hatch, 12/13/2013, hhatch.com, harold@hhatch.com
 *aa
 *
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
  WLTMMC mc(rstFileName.str().c_str());

//  // patch
//  mc.c()->collectInit(newfCollect);
//  AnalyzeMonkeyPatch patch(mc.space(), mc.pair());
//  mc.initAnalyze(&patch);
 
  // run sweeps
  mc.runNumSweepsRestart(100, rstFileName.str().c_str());
}
