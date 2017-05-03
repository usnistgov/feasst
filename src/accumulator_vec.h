/**
 * \file
 *
 * \brief
 *
 */

#ifndef ACCUMULATORVEC_H_
#define ACCUMULATORVEC_H_

#include "./accumulator.h"

namespace feasst {

class AccumulatorVec : public Accumulator {
 public:
  AccumulatorVec();
  virtual ~AccumulatorVec() {}

  /// accumulate values
  void accumulate(const int index, const double value);

  /// reset accumulator
  void reset();

  /// return normalized histogram of likelihood to accumulate a given index
  vector<double> hist() const;

  /// return total sum of all elements
  double sum() const;

  /// return the average of the average
  Accumulator average() const;
  double avOfAv() const { return average().average(); }

  /// read-only access of private data-members
  vector<Accumulator> vec() const { return accVec_; }
  Accumulator vec(const int i) const { return accVec_[i]; }
  int size() const { return static_cast<int>(accVec_.size()); }

 protected:
  vector<Accumulator> accVec_;    //<! vector of accumulators

  // error messaging
  void mout_(const char* messageType, std::ostream& message) {
    myOut(messageType, message, className_, verbose_);
  }
};

}  // namespace feasst

#endif  // ACCUMULATORVEC_H_

