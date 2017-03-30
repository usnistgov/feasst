/**
 * \file
 *
 * \brief interface for random number generator
 *
 */

#ifndef RANDOM_H_
#define RANDOM_H_

#include "base.h"

class Random : public Base {
 public:
  explicit Random(const unsigned long long iseed);
  explicit Random(const char* fileName);
  virtual ~Random() {}
  virtual void writeRestart(const char* fileName) = 0;

  /// seed random number generator
  virtual void seed(const unsigned long long iseed);

  /// uniform random doubleprecision number between 0 and 1
  virtual double uniform() = 0;

  /// uniform random integer between range min and max, inclusive
  virtual int uniform(const int min, const int max) = 0;

  /// random 32 bit unsigned integer
  unsigned int int32() { return (unsigned int) int64(); }

  /// random 64 bit integer
  virtual unsigned long long int64() = 0;

 protected:
  unsigned long long seed_;

  /// error messaging
  void mout_(const char* messageType, std::ostream& message) {
    myOut(messageType, message, className_, verbose_);
  }
};

#endif  // RANDOM_H_

