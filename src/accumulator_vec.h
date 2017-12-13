/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef ACCUMULATORVEC_H_
#define ACCUMULATORVEC_H_

#include "./accumulator.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Accumulate a series of values as a function of an order parameter.
 * In particular, this is the 1D-vector/array version of Accumulator.
 * HWH NOTE: This class is depreciated as it seems to have some duplication
 * with Histogram.
 */
class AccumulatorVec : public Accumulator {
 public:
  /// Constructor.
  AccumulatorVec();
  virtual ~AccumulatorVec() {}

  /// Add a value to the running sum of values and sum of squared values.
  /// @param [in] index vector index of accumulator
  void accumulate(const int index, const double value);

  /// Reset accumulator.
  void reset();

  /// Return normalized histogram of likelihood to accumulate a given index.
  vector<double> hist() const;

  /// Return total sum of all elements.
  double sum() const;

  /// Return the average of the average.
  Accumulator average() const;
  double avOfAv() const { return average().average(); }

  /// read-only access of private data-members
  vector<Accumulator> vec() const { return accVec_; }
  Accumulator vec(const int bin) const;
  int size() const { return static_cast<int>(accVec_.size()); }

 protected:
  vector<Accumulator> accVec_;    //<! vector of accumulators
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // ACCUMULATORVEC_H_

