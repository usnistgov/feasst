#include "feasst.h"

int main() {  // LJ, NVT_EXAMPLE
  feasst::Space space(3,    // 3D space
    {{"boxLength", "8"}});  // cubic PBC
  feasst::PairLJ pair(&space,
    {{"rCut", "3"},         // potential truncation
     {"cutType", "lrc"}});  // long range corrections
  feasst::CriteriaMetropolis criteria(1.2);  // beta = 1/k_B/T
  feasst::MC mc(&space, &pair, &criteria);
  feasst::transformTrial(&mc,
    {{"type", "translate"},
     {"maxMoveParam", "0.1"}});  // maximum displacement for each dimension
  mc.nMolSeek(50);  // add particles
  mc.initLog("log", 1e4);
  mc.initMovie("movie", 1e4);
  mc.runNumTrials(1e6);
}
