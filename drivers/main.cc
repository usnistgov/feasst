#include "pair_lj.h"
#include "mc.h"
#include "ui_abbreviated.h"

using namespace feasst;

int main() {
  Space space(3, 0);
  for (int dim = 0; dim < space.dimen(); ++dim) {
    space.lset(8, dim); // 8 box length
  }
  space.addMolInit("../../forcefield/data.lj");
  PairLJ pair(&space, 3);   // potential truncation at 3
  pair.initEnergy();
  CriteriaMetropolis criteria(1.2, 1.);  // 1/kT = 1.2
  MC mc(&space, &pair, &criteria);
  transformTrial(&mc, "translate", 0.1);
  mc.nMolSeek(50, "../../forcefield/data.lj");  // add 50 particles
  mc.initLog("log", 1e4);
  mc.initMovie("movie", 1e4);
  mc.runNumTrials(1e6);
}
