/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_field_slit_lj.h"

namespace feasst {

PairFieldSlitLJ::PairFieldSlitLJ(Space * space)
  : Pair(space),
    PairField(space),
    PairFieldSlit(space),
    PairFieldLJ(space) {
  defaultConstruction_();
}

PairFieldSlitLJ::PairFieldSlitLJ(Space * space,
  const char* fileName)
  : Pair(space, fileName),
    PairField(space, fileName),
    PairFieldSlit(space, fileName),
    PairFieldLJ(space, fileName) {
  defaultConstruction_();
}

void PairFieldSlitLJ::defaultConstruction_() {
  className_.assign("PairFieldSlitLJ");
}

void PairFieldSlitLJ::writeRestart(const char* fileName) {
  PairFieldSlit::writeRestart(fileName);
  PairFieldLJ::writeRestart(fileName);
}

void PairFieldSlitLJ::fieldSiteSlit_(const int &siteType, double * energy,
  double * force, const double &dist) {
  if (dist <= 0.) {
    *energy = NUM_INF;
  } else {
    *energy = eps_[siteType]*pow(sig_[siteType]/(dist + delta_), alpha_);
  }
}

shared_ptr<PairFieldSlitLJ> makePairFieldSlitLJ(std::shared_ptr<Space> space) {
  return make_shared<PairFieldSlitLJ>(space.get());
}

}  // namespace feasst
