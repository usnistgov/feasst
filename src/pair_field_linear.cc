/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_field_linear.h"

namespace feasst {

PairFieldLinear::PairFieldLinear(Space * space, const char* fileName)
  : Pair(space, fileName) {
  // makeFactory requires space argument in Pair classes
  if (space == NULL) {} // remove unused parameter warning
  defaultConstruction_();
}

void PairFieldLinear::defaultConstruction_() {
}

void PairFieldLinear::writeRestart(const char* fileName) {
  std::ofstream file(fileName, std::ios_base::app);
}

double PairFieldLinear::linear_energy_(const int &siteType, const double &dist) {
  const double halfSig = 0.5*sig_[siteType];
  const double rCut = rCutij_[siteType][siteType];
  if (dist < halfSig) {
    return NUM_INF;
  } else if (dist <= rCut) {
    return eps_[siteType]*(dist/rCut - 1.);
  }
  return 0.;
}

}  // namespace feasst
