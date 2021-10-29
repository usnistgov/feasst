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
#include <string>
#include "utils/include/arguments.h"

namespace feasst {

/**
 * Accumulate a series of values to compute the average, standard deviation,
 * and the standard deviation of the average of uncorrelated data using the
 * blocking method.
 * Flyvbjerg, Error estimates on averages of correlated data,
 * http://dx.doi.org/10.1063/1.457480
 * Also, see Frenkel and Smit, Understanding Molecular Simulations, Appendix D
 * and Case study 4 (pages 98-100, Figure 4.4).
 */
class Accumulator {
 public:
  /**
    args:
    - num_moments: maximum number of moments (default: 5).
    - max_block_operations: maximum number of blocking operations (default: 6).
   */
  explicit Accumulator(argtype args = argtype());
  explicit Accumulator(argtype * args);

  /// Return the maximum number of block operations.
  int max_block_operations() const { return max_block_operations_; }

  /// Return the maximum number of moments.
  int num_moments() const { return static_cast<int>(val_moment_.size()); }

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
    Return standard deviation of the block average.
    Return 0 if not enough data.
   */
  double block_stdev(
    /// Use the given operation.
    /// If -1, have the maximum stdev over all blocks.
    const int num_op = -1,
    /// If num_op == -1, minimum number of blocks to consider stdev.
    const int min_blocks = 10) const;

  /// Return standard deviation of the standard deviation of the block average.
  double block_std_of_std(const int num_op = 0) const;

  /// Return the block averages
  const std::vector<std::shared_ptr<Accumulator> >& block_averages() const {
    return block_averages_; }

  /// Return the size of each block.
  const std::vector<double>& block_size() const { return block_size_; }

  /// Return number of values accumulated.
  long double num_values() const { return val_moment_[0]; }

  /// Return sum of all values accumulated.
  long double sum() const { return val_moment_[1]; }

  /// Same as above, but truncated to double precision for Python interface.
  double sum_dble() const { return static_cast<double>(sum()); }

  /// Return sum of the square of all values accumulated.
  long double sum_of_squared() const { return val_moment_[2]; }

  /// Same as above, but truncated to double precision for Python interface.
  double sum_of_squared_dble() const { return static_cast<double>(val_moment_[2]); }

  /// Zero all accumulated values.
  void reset();

  /// Return the maximum value accumulated.
  double max() const { return max_; }

  /// Return the minimum value accumulated.
  double min() const { return min_; }

  /// Return the moments as a double.
  /// This is the sum of the value^index.
  /// Thus, moment(0) is the number of accumulated values.
  /// And moment(1) is the sum of the values.
  /// moment(2) is the sum of the squared values, etc.
  double moment(const int index) const {
    return static_cast<double>(val_moment_[index]); }

  /// Return the moments.
  std::vector<long double> moments() const { return val_moment_; }

  // HWH: requires further testing/documentation
  // moment = 1/N sum((i - av)^n)
  // Return the n-th central moment.
  // https://en.wikipedia.org/wiki/Central_moment
  double central_moment(const int n);

  /// Return the header of the human readable status.
  std::string status_header() const;

  /// Return human readable status.
  std::string status() const;

  /// Combine header and status.
  std::string str() const;

  /// Return the last value accumulated.
  double last_value() const;

  /// Return true if the Accumulator is equivalent within the given
  /// confidence interval based on the block average standard deviations.
  bool is_equivalent(const Accumulator& accumulator,
    const double t_factor,
    const int num_op = 0,  // number of block operations
    const bool verbose = false) const;

  void serialize(std::ostream& ostr) const;
  explicit Accumulator(std::istream& istr);

 private:
  double max_, min_;
  double last_value_;
  std::vector<long double> val_moment_;

  // block averaging variables
  int block_power_ = 2;
  int max_block_operations_;
  std::vector<double> block_size_;
  std::vector<long double> sum_block_;
  std::vector<std::shared_ptr<Accumulator> > block_averages_;

  // Set the highest order of moments recorded.
  void set_moments_(const int num_moments);
};

inline std::shared_ptr<Accumulator> MakeAccumulator(argtype args = argtype()) {
  return std::make_shared<Accumulator>(args); }

}  // namespace feasst

#endif  // FEASST_MATH_ACCUMULATOR_H_
