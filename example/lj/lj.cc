/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include "pair_lj.h"
#include "mc.h"
#include "trial_transform.h"

int main() {
  feasst::Space space(3, 0);
  space.lset(8);  // 8 box length cubic PBC
  space.addMolInit("../../forcefield/data.lj");
  feasst::PairLJ pair(&space, 3);   // potential truncation at 3
  pair.initEnergy();
  feasst::CriteriaMetropolis criteria(1.2, 1.);  // 1/kT = 1.2
  feasst::MC mc(&space, &pair, &criteria);
  feasst::transformTrial(&mc, "translate", 0.1);
  mc.nMolSeek(50, "../../forcefield/data.lj");  // add 50 particles
  mc.initLog("log", 1e4);
  mc.initMovie("movie", 1e4);
  mc.runNumTrials(1e6);
}
