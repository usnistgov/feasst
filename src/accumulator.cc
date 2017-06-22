#include "./accumulator.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

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

void Accumulator::accumulate(double value) {
  ++nValues_;
  sum_ += value;
  sumSq_ += value*value;

  if (max_ < value) max_ = value;
  if (min_ > value) min_ = value;

  // accumulate block averages
  if (nBlock_ != 0) {
    sumBlock_ += value;
    if (nValues_ % nBlock_ == 0) {
      if (blockAvs_ == NULL) blockAvs_ = make_shared<Accumulator>();
      blockAvs_->accumulate(sumBlock_ / static_cast<double>(nBlock_));
      sumBlock_ = 0.;
    }
  }
}

void Accumulator::reset() {
  sum_ = 0;
  nValues_ = 0;
  sumSq_ = 0;
  setBlock();
  sumBlock_ = 0;
  max_ = -NUM_INF;
  min_ = NUM_INF;
}

double Accumulator::average() const {
  double av = 0;
  if (nValues_ > 0) {
    av = static_cast<double>(sum_ / static_cast<double>(nValues_));
  }
  return av;
}

double Accumulator::stdev() const {
  double stdev = 0;
  if (nValues_ > 1) {
    const double fluct = sumSq_/static_cast<double>(nValues_)
                       - pow(average(), 2);
    if (fluct > 0.) {
      stdev = sqrt(fluct*nValues_ / static_cast<double>(nValues_ - 1));
    }
  }
  return stdev;
}

double Accumulator::blockStdev() const {
  if (nBlock_ != 0) {
    if (blockAvs_ != NULL) {
      if (blockAvs_->nValues() > 1) {
        return blockAvs_->stdev()/sqrt(static_cast<int>(blockAvs_->nValues()));
      }
    }
  }
  return 0;
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

