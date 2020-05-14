/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_field_lj.h"

namespace feasst {

PairFieldLJ::PairFieldLJ(Space * space, const char* fileName)
 : Pair(space, fileName) {
  // makeFactory requires space argument in Pair classes
  if (space == NULL) {} // remove unused parameter warning
  defaultConstruction_();
  alpha_ = fstoi("confineAlpha", fileName);
  delta_ = fstod("confineDelta", fileName);
}

void PairFieldLJ::defaultConstruction_() {
  initAlpha();
  initDelta();
}

void PairFieldLJ::writeRestart(const char* fileName) {
  std::ofstream file(fileName, std::ios_base::app);
  file << "# confineAlpha " << alpha_ << endl;
  file << "# confineDelta " << delta_ << endl;
}

}  // namespace feasst
