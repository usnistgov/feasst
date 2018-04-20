/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef ACCUMULATOR_H_
#define ACCUMULATOR_H_

#include <vector>
#include <memory>
#include "./base.h"

namespace feasst {

/**
 * Accumulate a series of values to compute the average, standard deviation,
 * and the standard deviation of the average of uncorrelated data using the
 * blocking method.
 */
class Accumulator : public Base {
 public:
  /// Constructor
  Accumulator();

  /// Add a value to the running sum of values and higher moments.
  virtual void accumulate(double value);

  // Dummy virtual function to avoid warnings for hidden overloaded virtual.
  virtual void accumulate(const int index, const double value) {
    ASSERT(0, "This function should not be used");
  }

  /// Return average of accumulated values.
  double average() const;

  /// Return standard deviation e.g., fluctuation of all (correlated) values.
  double std() const;

  /// Return the standard deviation of the average (e.g., std/sqrt(num_samples))
  double stdOfAv() const { return std()/sqrt(static_cast<double>(nValues())); }

  /**
   * Set the size of a block to compute (uncorrelated) standard errors of the
   * block averages.
   *
   * Flyvbjerg, Error estimates on averages of correlated data,
   * http://dx.doi.org/10.1063/1.457480
   */
  void setBlock(
    const long long nBlock = 1e5  //!< nBlock number of values per block
    ) { nBlock_ = nBlock; }

  /// Return standard deviation of the block averages (0 if not enough blocks).
  double blockStdev() const;

  /// Constructor for checkpointing. Note that block averages are not saved.
  Accumulator(const long long nValues, const long double sum,
              const long double sumSq);

  /// Return number of values accumulated.
  long long nValues() const { return (long long)nValues_; }

  /// Return sum of all values accumulated.
  long double sum() const { return sum_; }

  /// Return sum of the squared of all values accumulated.
  double sumDble() const { return sum_; }

  /// Return sum of the square of all values accumulated.
  long double sumSq() const { return sumSq_; }

  /// Return sum of the squared of all values accumulated as double.
  double sumSqDble() const { return sumSq_; }

  /// Zero all accumulated values.
  virtual void reset();

  /// Return number of values per block.
  long long nBlock() const { return nBlock_; }

  /// Return the maximum value accumulated.
  double max() const { return max_; }

  /// Return the minimum value accumulated.
  double min() const { return min_; }

  /// Set the highest order of moments recorded.
  void initMoments(const int nMoments = 2) { valMoment_.resize(nMoments); }

  /// Return the moments. Note that this is not supported with checkpointing.
  vector<long double> valMoment() const { return valMoment_; }

  virtual ~Accumulator() {}

 protected:
  long long nValues_;
  long double sum_;
  long double sumSq_;
  double max_, min_;

  /// maximum number of moments before Taylor series truncation
  ///  e.g., 2 (default) includes moments 1 and 2
  vector<long double> valMoment_;

  // block averaging variables
  long long nBlock_;
  long double sumBlock_;         //!< sum of all values accumulated
  shared_ptr<Accumulator> blockAvs_;    //!< accumulate averages of each block
};

}  // namespace feasst

#endif  // ACCUMULATOR_H_

