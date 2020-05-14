/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_field_slit.h"

namespace feasst {

PairFieldSlit::PairFieldSlit(Space * space)
  : Pair(space),
    PairField(space) {
  defaultConstruction_();
}

PairFieldSlit::PairFieldSlit(Space * space, const char* fileName)
  : Pair(space, fileName),
    PairField(space, fileName) {
  defaultConstruction_();
  confine_dimen_ = fstoi("confineDimen", fileName);
  upper_ = fstoi("slitUpper", fileName);
  lower_ = fstoi("slitLower", fileName);
}

void PairFieldSlit::defaultConstruction_() {
  className_.assign("PairFieldSlit");
  initSlit();
}

void PairFieldSlit::writeRestart(const char* fileName) {
  PairField::writeRestart(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# confineDimen " << confine_dimen_ << endl;
  file << "# slitUpper " << upper_ << endl;
  file << "# slitLower " << lower_ << endl;
}

void PairFieldSlit::fieldSite_(const int &siteType, double * energy, double * force,
  const double &x, const double &y, const double &z) {
  * energy = 0;
  double pos = NUM_INF;
  if (confine_dimen_ == 0) {
    pos = x;
  } else if (confine_dimen_ == 1) {
    pos = y;
  } else if (confine_dimen_ == 2) {
    pos = z;
  }

  // compute energy of the upper wall
  double energyPerWall;
  fieldSiteSlit_(siteType, &energyPerWall, force, upper_ - pos);
  *energy += energyPerWall;

  // compute energy of the lower wall
  fieldSiteSlit_(siteType, &energyPerWall, force, pos - lower_);
  *energy += energyPerWall;
}

void PairFieldSlit::initSlit(const int dimen, const double upper,
  const double lower) {
  confine_dimen_ = dimen;
  upper_ = upper;
  lower_ = lower;
}

shared_ptr<PairFieldSlit> makePairFieldSlit(std::shared_ptr<Space> space) {
  return make_shared<PairFieldSlit>(space.get());
}

}  // namespace feasst
