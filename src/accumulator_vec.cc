/**
 * \file
 *
 * \brief
 */

#include "./accumulator_vec.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Constructor
 */
AccumulatorVec::AccumulatorVec() {
  verbose_ = 0;
  className_.assign("AccumulatorVec");
  reset();
}

/**
 * accumulate values
 */
void AccumulatorVec::accumulate(
  const int index,  //!< vector index of accumulator
  const double value) {  //!< value to accumulate
  // resize vector if index is out of range
  if (static_cast<int>(accVec_.size()) <= index) accVec_.resize(index+1);
  accVec_[index].accumulate(value);
}

/**
 * reset accumulator
 */
void AccumulatorVec::reset() {
  accVec_.clear();
}

/**
 * return normalized histogram of likelihood to accumulate a given index
 */
vector<double> AccumulatorVec::hist() const {
  long long n = 0;
  for (int i = 0; i < size(); ++i) {
    n += accVec_[i].nValues();
  }
  vector<double> p(size());
  for (int i = 0; i < size(); ++i) {
    p[i] = accVec_[i].nValues() / static_cast<double>(n);
  }
  return p;
}

/**
 * return total sum of all elements
 */
double AccumulatorVec::sum() const {
  long double sum = 0;
  for (unsigned int i = 0; i < accVec_.size(); ++i) {
    sum += accVec_[i].sum();
  }
  return sum;
}

/// return the average of the average
Accumulator AccumulatorVec::average() const {
  Accumulator Acc;
  for (unsigned int i = 0; i < accVec_.size(); ++i) {
    Acc.accumulate(accVec_[i].average());
  }
  return Acc;
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_






