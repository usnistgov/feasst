/**
 * \file
 *
 * \brief
 *
 */

#ifndef ACCUMULATOR_H_
#define ACCUMULATOR_H_

#include "./base.h"

class Accumulator : public Base {
 public:
  Accumulator();
  Accumulator(const long long nValues, const long double sum,
              const long double sumSq);
  virtual ~Accumulator() {}

  /// accumulate values
  virtual void accumulate(double i);

  /// reset accumulator
  virtual void reset();

  /// running average
  double average() const;

  /// running standard deviation
  double stdev() const;

  /// set block size
  void setBlock(const long long nBlock) { nBlock_ = nBlock; }

  /// output block average stdev
  double blockStdev() const;

  /// read-only access of private data-members
  long long nValues() const { return (long long)nValues_; }
  long double sum() const { return sum_; }
  long double sumSq() const { return sumSq_; }

 protected:
  long long nValues_;   //!< number of values accumulated
  long double sum_;         //!< sum of all values accumulated
  long double sumSq_;       //!< squared sum of all values accumulated

  // block averaging variables
  long long nBlock_;    //!< number of values per block
  long double sumBlock_;         //!< sum of all values accumulated
  shared_ptr<Accumulator> blockAvs_;    //!< accumulate averages of each block

  // error messaging
  void mout_(const char* messageType, std::ostream& message) {
    myOut(messageType, message, className_, verbose_);
  }
};

#endif  // ACCUMULATOR_H_

