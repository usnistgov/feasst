/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef FEASST_MATH_ACCUMULATOR_H_
#define FEASST_MATH_ACCUMULATOR_H_

#include <vector>
#include <memory>
#include "utils/include/utils_io.h"

namespace feasst {

/**
 * Accumulate a series of values to compute the average, standard deviation,
 * and the standard deviation of the average of uncorrelated data using the
 * blocking method.
 */
class Accumulator {
 public:
  Accumulator();

  /// Add a value to the running sum of values and higher moments.
  void accumulate(double value);

  /// Return average of accumulated values.
  double average() const;

  /// Return standard deviation e.g., fluctuation of all (correlated) values.
  double stdev() const;

  /// Same as above.
  double std() const { return stdev(); }

  /// Return the standard deviation of the average (e.g., std/sqrt(num_samples))
  double stdev_of_av() const;

  /**
   * Set the size of a block to compute (uncorrelated) standard errors of the
   * block averages.
   *
   * Flyvbjerg, Error estimates on averages of correlated data,
   * http://dx.doi.org/10.1063/1.457480
   */
  void set_block(
    /// number of values per block
    const long long num_block = 1e5);

  /// Return standard deviation of the block averages (0 if not enough blocks).
  double block_stdev() const;

  /// Return number of values accumulated.
  long long num_values() const { return (long long)num_values_; }

  /// Return sum of all values accumulated.
  long double sum() const { return sum_; }

  /// Return sum of the square of all values accumulated.
  long double sum_of_squared() const { return sum_squared_; }

  /// Zero all accumulated values.
  void reset();

  /// Return number of values per block.
  long long num_block() const { return num_blocks_; }

  /// Return the maximum value accumulated.
  double max() const { return max_; }

  /// Return the minimum value accumulated.
  double min() const { return min_; }

  /// Set the highest order of moments recorded.
  void set_moments(const int num_moments = 2);

  /// Return the moments.
  std::vector<long double> moments() const { return val_moment_; }

  /// Return human readable status as a string.
  std::string str() const;

  /// Return the last value accumulated.
  double last_value() const;

  void serialize(std::ostream& ostr) const;
  Accumulator(std::istream& istr);

 private:
  long long num_values_;
  long double sum_;
  long double sum_squared_;
  double max_, min_;
  double last_value_;

  /// maximum number of moments before Taylor series truncation
  ///  e.g., 2 (default) includes moments 1 and 2
  std::vector<long double> val_moment_;

  // block averaging variables
  long long num_blocks_;
  long double sum_block_;         //!< sum of all values accumulated

  /// accumulate averages of each block
  std::shared_ptr<Accumulator> block_averages_;
};

}  // namespace feasst

#endif  // FEASST_MATH_ACCUMULATOR_H_
