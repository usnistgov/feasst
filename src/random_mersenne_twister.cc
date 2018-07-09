/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "random_mersenne_twister.h"

namespace feasst {

RandomMersenneTwister::RandomMersenneTwister(const unsigned long long iseed
  ) : Random(iseed) {
  defaultConstruction_();
  seed(iseed);
}

RandomMersenneTwister::RandomMersenneTwister(const char* fileName
  ) : Random(fileName) {
  defaultConstruction_();
  seed_ = fstoull("seed", fileName);
}

void RandomMersenneTwister::defaultConstruction_() {
  verbose_ = 0;
  className_.assign("RandomMersenneTwister");
  dis_double_ = std::uniform_real_distribution<double>(0.0, 1.0);
}

void RandomMersenneTwister::seed(const unsigned long long iseed) {
  Random::seed(iseed);
  generator_ = std::mt19937(iseed);
}

double RandomMersenneTwister::uniform() {
  return dis_double_(generator_);
}

unsigned long long RandomMersenneTwister::int64() {
  return dis_longlong_(generator_);
}

void RandomMersenneTwister::writeRestart(const char* fileName) {
  fileBackUp(fileName);
  std::ofstream file(fileName);
  file << "# seed " << seed_ << endl;
}

}  // namespace feasst
