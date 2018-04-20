/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef RANDOM_NR3_H_
#define RANDOM_NR3_H_

#include "random.h"

namespace feasst {

/**
 * Numerical recipes third edition random number generator.
 */
class RandomNR3 : public Random {
 public:
  RandomNR3(const unsigned long long seed);
  RandomNR3(const char* fileName);

  // Overloaders for virtual functions. See base class for comments.
  ~RandomNR3() {};
  void writeRestart(const char* fileName);
  void seed(const unsigned long long seed);
  double uniform();
  unsigned long long int64();

 protected:
  /// NR3 specific stored state variables u, v and w.
  unsigned long long u_, v_, w_;
};

}  // namespace feasst

#endif  // RANDOM_NR3_H_

