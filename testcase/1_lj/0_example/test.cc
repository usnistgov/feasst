#include "feasst.h"

int main() {  // LJ, NVT_EXAMPLE
  feasst::Space space(3);   // 3D space
  space.initBoxLength(8);   // cubic PBC
  feasst::PairLJ pair(&space, 3,  // potential truncation
    {{"cutType", "lrc"}});
  feasst::CriteriaMetropolis criteria(1.2,  // beta = 1/k_B/T
    1.);  // lnz = beta mu - 3ln(Lambda)
  feasst::MC mc(&space, &pair, &criteria);
  const double maxMoveParam = 0.1;  // maximum displacement for each dimension
  feasst::transformTrial(&mc, "translate", maxMoveParam);
  mc.nMolSeek(50);  // add particles
  mc.initLog("log", 1e4);
  mc.initMovie("movie", 1e4);
  mc.runNumTrials(1e6);
}
