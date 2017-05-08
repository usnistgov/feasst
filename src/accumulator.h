/**
 * This class is used to accumulate averages and standard deviations.
 */

#ifndef ACCUMULATOR_H_
#define ACCUMULATOR_H_

#include "./base.h"

namespace feasst {

class Accumulator : public Base {
 public:
  Accumulator();    //!< Constructor

  /**
   * Accumulator constructor for checkpointing.
   */
  Accumulator(const long long nValues, const long double sum,
              const long double sumSq);

  virtual ~Accumulator() {}

  /// accumulate values
  virtual void accumulate(double value);

  /// reset accumulator
  virtual void reset();

  /// \return average of accumulated values
  double average() const;

  /// \return standard deviation e.g., fluctuation of all (correlated) values
  double stdev() const;

  /*!
   * set the size of a block to compute (uncorrelated) standard errors of the
   * block averages
   *
   * Flyvbjerg, Error estimates on averages of correlated data,
   * http://dx.doi.org/10.1063/1.457480
   *
   * \param nBlock number of values per block
   */
  void setBlock(const long long nBlock = 1e5) { nBlock_ = nBlock; }

  /// \return standard deviation of the block averages (0 if not enough blocks)
  double blockStdev() const;

  /// number of values accumulated
  long long nValues() const { return (long long)nValues_; }

  /// sum of all values accumulated
  long double sum() const { return sum_; }

  /// sum of the square of all values accumulated
  long double sumSq() const { return sumSq_; }

  /// number of values per block
  long long nBlock() const { return nBlock_; }

 protected:
  long long nValues_;
  long double sum_;
  long double sumSq_;

  // block averaging variables
  long long nBlock_;
  long double sumBlock_;         //!< sum of all values accumulated
  shared_ptr<Accumulator> blockAvs_;    //!< accumulate averages of each block
};

}  // namespace feasst

#endif  // ACCUMULATOR_H_

