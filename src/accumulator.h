#ifndef ACCUMULATOR_H_
#define ACCUMULATOR_H_

#include "./base.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Accumulate a series of values to compute the average, standard deviation,
 * and the standard deviation of the average of uncorrelated data using the
 * blocking method.
 */
class Accumulator : public Base {
 public:
  Accumulator();
  virtual ~Accumulator() {}

  /// Constructor for checkpointing. Note that block averages are not saved.
  Accumulator(const long long nValues, const long double sum,
              const long double sumSq);

  /// Return number of values accumulated.
  long long nValues() const { return (long long)nValues_; }

  /// Sum of all values accumulated.
  long double sum() const { return sum_; }
  double sumDble() const { return sum_; }

  /// Sum of the square of all values accumulated.
  long double sumSq() const { return sumSq_; }
  double sumSqDble() const { return sumSq_; }

  /// Add a value to the running sum of values and sum of squared values.
  virtual void accumulate(double value);

  /// Zero all accumulated values.
  virtual void reset();

  /// Return average of accumulated values.
  double average() const;

  /// Return standard deviation e.g., fluctuation of all (correlated) values.
  double stdev() const;

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

  /// Return number of values per block.
  long long nBlock() const { return nBlock_; }

  /// Return standard deviation of the block averages (0 if not enough blocks).
  double blockStdev() const;

 protected:
  long long nValues_;
  long double sum_;
  long double sumSq_;

  // block averaging variables
  long long nBlock_;
  long double sumBlock_;         //!< sum of all values accumulated
  shared_ptr<Accumulator> blockAvs_;    //!< accumulate averages of each block
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // ACCUMULATOR_H_

