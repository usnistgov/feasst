#ifndef RANDOM_NR3_H_
#define RANDOM_NR3_H_

#include "random.h"

namespace feasst {

/**
 * Numerical recipes third edition random number generator.
 */
class RandomNR3 : public Random {
 public:
  /// Construct from seed.
  RandomNR3(const unsigned long long seed);

  /// Construct from checkpoint file.
  RandomNR3(const char* fileName);

  virtual ~RandomNR3() {};

  /// Write checkpoint file.
  virtual void writeRestart(const char* fileName);

  /// Seed random number generator.
  void seed(const unsigned long long seed);

  /// Return uniform random number between 0 and 1.
  double uniform();

  // NOTE HWH: Depreciate this
  // uniform random integer between range min and max, inclusive
  int uniform(const int min, const int max);

  /// Return random 64 bit integer.
  unsigned long long int64();

 protected:
  unsigned long long u_, v_, w_;
};

}  // namespace feasst

#endif  // RANDOM_NR3_H_

