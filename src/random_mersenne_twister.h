/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef RANDOM_MERSENNE_TWISTER_H_
#define RANDOM_MERSENNE_TWISTER_H_

#include <random>
#include "random.h"

namespace feasst {

/**
 * STL C++ implementation of the Mersenne Twister 19937 random number generator.
 */
class RandomMersenneTwister : public Random {
 public:
  RandomMersenneTwister(const unsigned long long seed);
  RandomMersenneTwister(const char* fileName);

  // Overloaders for virtual functions. See base class for comments.
  ~RandomMersenneTwister() {};
  void writeRestart(const char* fileName);
  void seed(const unsigned long long seed);
  double uniform();
  unsigned long long int64();

 protected:
  std::uniform_real_distribution<double> dis_double_;
  std::uniform_int_distribution<unsigned long long> dis_longlong_;
  std::mt19937 generator_;
  void defaultConstruction_();
};

}  // namespace feasst

#endif  // RANDOM_MERSENNE_TWISTER_H_
