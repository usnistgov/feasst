/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./accumulator.h"

namespace feasst {

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
  valMoment_[0] = sum_;
  valMoment_[1] = sumSq_;
}

void Accumulator::accumulate(double value) {
  ++nValues_;
  sum_ += value;
  sumSq_ += value*value;

  // accumulate moments
  double valmo = 1.;
  for (int mo = 0; mo < static_cast<int>(valMoment_.size()); ++mo) {
    valmo *= value;
    valMoment_[mo] += valmo;
  }

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
  initMoments();
}

double Accumulator::average() const {
  double av = 0;
  if (nValues_ > 0) {
    av = static_cast<double>(sum_ / static_cast<double>(nValues_));
  }
  return av;
}

double Accumulator::std() const {
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
        return blockAvs_->std()/sqrt(static_cast<int>(blockAvs_->nValues()));
      }
    }
  }
  return 0;
}

}  // namespace feasst

