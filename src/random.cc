/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "random.h"

namespace feasst {

Random::Random(const unsigned long long iseed
  ) : seed_(iseed) {
  verbose_ = 0;
  className_.assign("Random");
}

Random::Random(const char* fileName) {
  verbose_ = 0;
  className_.assign("Random");
  (void) fileName;  // avoid unused parameter warning
}

void Random::seed(const unsigned long long iseed) {
  seed_ = iseed;
}

}  // namespace feasst



