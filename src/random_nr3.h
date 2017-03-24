/**
 * \file
 *
 * \brief numerical recipes third edition random number generator
 *
 */

#ifndef RANDOMNR3_H_
#define RANDOMNR3_H_

#include "random.h"

class RandomNR3 : public Random {
 public:
  RandomNR3(const unsigned long long seed);
  RandomNR3(const char* fileName);
  virtual ~RandomNR3() {};
  virtual void writeRestart(const char* fileName);

  // seed random number generator
  void seed(const unsigned long long seed);

  // uniform random doubleprecision number between 0 and 1
  double uniform();

  // uniform random integer between range min and max, inclusive
  int uniform(const int min, const int max);

  // random 64 bit integer
  unsigned long long int64();

 protected:
  unsigned long long u_, v_, w_;
};

#endif  // RANDOMNR3_H_

