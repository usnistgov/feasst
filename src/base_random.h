/**
 * \file
 *
 * \brief interface for all classes to inherit
 *
 */

#ifndef BASE_RANDOM_H_
#define BASE_RANDOM_H_

#include "./base.h"

class Random;

class BaseRandom : public Base {
 public:
  BaseRandom();
  explicit BaseRandom(const char* fileName);
  virtual ~BaseRandom();

  /// write restart file for random number generator, if initialized
  void writeRngRestart(const char* fileName);

  /// derived objects may preform additional reconstruction
  void reconstructDerived();

  /// initialize random number generator
  void initRNG(Random *ran);
  void initRNG(unsigned long long seed);
  void initRNG();
  void initRNG(const char* fileName);

  /// clear random number generator
  void clearRNG();

  /// uniform random doubleprecision number between 0 and 1
  double uniformRanNum();

  /// standard normal random doubleprecision number
  double stdNormRanNum();

  /// multivariate gaussian random number of given variance and dimension
  vector<double> stdRanNum(const double variance, const int dimen);

  /// gaussian distribution random number
  double gaussRanNum(const double variance, const double av) {
    return av + variance*stdNormRanNum();
  }

  /// uniform random integer between range min and max, inclusive
  int uniformRanNum(const int min, const int max);

  /// gassuian distribution random number from Frenkel and Smit
  double gaussRanNumFS();
 
 protected:
  Random* ranNum_;      //!< pointer to random number generator

  /// set true if the class owns the pointer to the random number generator
  bool ranNumOwned_;
};

#endif  // BASE_RANDOM_H_

