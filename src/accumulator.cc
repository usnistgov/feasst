/**
 * \file
 *
 * \brief
 */

#include "./accumulator.h"

/**
 * Constructor
 */
Accumulator::Accumulator() {
  verbose_ = 0;
  className_.assign("Accumulator");
  reset();
}
Accumulator::Accumulator(const long long nValues, const long double sum,
  const long double sumSq) {
  reset();
  nValues_ = nValues;
  sum_ = sum;
  sumSq_ = sumSq;
}

/**
 * accumulate values
 */
void Accumulator::accumulate(double i) {
  ++nValues_;
  sum_ += i;
  sumSq_ += i*i;

  // accumulate block averages
  if (nBlock_ != 0) {
    sumBlock_ += i;
    if (nValues_ % nBlock_ == 0) {
      if (blockAvs_ == NULL) blockAvs_ = make_shared<Accumulator>();
      blockAvs_->accumulate(sumBlock_ / static_cast<double>(nBlock_));
      sumBlock_ = 0.;
    }
  }
}

/**
 * reset accumulator
 */
void Accumulator::reset() {
  sum_ = 0;
  nValues_ = 0;
  sumSq_ = 0;
  nBlock_ = 1e5;
  sumBlock_ = 0;
}

/**
 * running average
 */
double Accumulator::average() const {
  double av = 0;
  if (nValues_ > 0) {
    av = static_cast<double>(sum_ / static_cast<double>(nValues_));
  }
  return av;
}

/**
 * running standard deviation
 *  this is not the standard error, but a measure of the fluctuation of the
 *  value about the average
 */
double Accumulator::stdev() const {
  double stdev = 0;
  if (nValues_ > 1) {
    stdev = sqrt((sumSq_/static_cast<double>(nValues_) - pow(average(), 2))
      * nValues_ / static_cast<double>(nValues_ - 1));
  }
  return stdev;
}

/**
 * output block average stdev
 */
double Accumulator::blockStdev() const {
  // cout << "hi " << nBlock_ << " " << blockAvs_ == NULL << endl;
  if (nBlock_ != 0) {
    if (blockAvs_ != NULL) {
      if (blockAvs_->nValues() > 1) {
        // cout << blockAvs_->nValues() << " " << blockAvs_->stdev() << endl;
        return blockAvs_->stdev()/sqrt(static_cast<int>(blockAvs_->nValues()));
      }
    }
  }
  return 0;
}
