#ifndef RANDOM_H_
#define RANDOM_H_

#include "base.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * General interface for random number generators.
 */
class Random : public Base {
 public:
  /// Construct from seed.
  explicit Random(const unsigned long long iseed);

  /// Construct from checkpoint file.
  explicit Random(const char* fileName);

  virtual ~Random() {}

  /// Write checkpoint file.
  virtual void writeRestart(const char* fileName) = 0;

  /// Seed random number generator.
  virtual void seed(const unsigned long long iseed);

  /// Return uniform random number between 0 and 1.
  virtual double uniform() = 0;

  /// Return uniform random integer between range min and max, inclusive.
  int uniform(const int min, const int max) {
    return int64() % (max - min + 1) + min;
    }; 

  // NOTE HWH: Depreciate this
  /// Return random 32 bit unsigned integer.
  unsigned int int32() { return (unsigned int) int64(); }

  /// Return random 64 bit integer.
  virtual unsigned long long int64() = 0;

 protected:
  unsigned long long seed_;
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // RANDOM_H_

