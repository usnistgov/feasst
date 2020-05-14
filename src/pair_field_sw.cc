/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_field_sw.h"

namespace feasst {

PairFieldSW::PairFieldSW(Space * space, const char* fileName)
  : Pair(space, fileName) {
  // makeFactory requires space argument in Pair classes
  if (space == NULL) {} // remove unused parameter warning
  defaultConstruction_();
}

void PairFieldSW::defaultConstruction_() {
}

void PairFieldSW::writeRestart(const char* fileName) {
  std::ofstream file(fileName, std::ios_base::app);
}

double PairFieldSW::sw_energy_(const int &siteType, const double &dist) {
  const double halfSig = 0.5*sig_[siteType];
  if (dist < halfSig) {
    return NUM_INF;
  } else if (dist <= rCutij_[siteType][siteType]) {
    return -eps_[siteType];
  }
  return 0.;
}

}  // namespace feasst
