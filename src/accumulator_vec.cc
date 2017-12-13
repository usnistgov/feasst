/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./accumulator_vec.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

AccumulatorVec::AccumulatorVec() {
  verbose_ = 0;
  className_.assign("AccumulatorVec");
  reset();
}

void AccumulatorVec::accumulate(
  const int index,
  const double value) {
  // resize vector if index is out of range
  if (static_cast<int>(accVec_.size()) <= index) accVec_.resize(index+1);
  accVec_[index].accumulate(value);
}

void AccumulatorVec::reset() {
  accVec_.clear();
}

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

double AccumulatorVec::sum() const {
  long double sum = 0;
  for (unsigned int i = 0; i < accVec_.size(); ++i) {
    sum += accVec_[i].sum();
  }
  return sum;
}

Accumulator AccumulatorVec::average() const {
  Accumulator Acc;
  for (unsigned int i = 0; i < accVec_.size(); ++i) {
    Acc.accumulate(accVec_[i].average());
  }
  return Acc;
}

Accumulator AccumulatorVec::vec(const int bin) const {
  if (accVec_.size() < 1) {
    Accumulator null;
    return null;
  }
  return accVec_[bin];
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_






