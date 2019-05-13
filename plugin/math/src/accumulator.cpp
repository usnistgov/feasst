/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

// #include <cmath>
#include <math.h>
#include <numeric>
#include "math/include/accumulator.h"
#include "math/include/constants.h"

namespace feasst {

Accumulator::Accumulator() {
  reset();
}

//Accumulator::Accumulator(const long long num_values, const long double sum,
//  const long double sum_squared) {
//  reset();
//  num_values_ = num_values;
//  sum_ = sum;
//  sum_squared_ = sum_squared;
//}

void Accumulator::accumulate(double value) {
  ++num_values_;
  sum_ += value;
  sum_squared_ += value*value;

  // accumulate moments
  double valmo = 1.;
  for (int mo = 0; mo < static_cast<int>(val_moment_.size()); ++mo) {
    valmo *= value;
    val_moment_[mo] += valmo;
  }

  if (max_ < value) max_ = value;
  if (min_ > value) min_ = value;

  // accumulate block averages
  if (num_blocks_ != 0) {
    sum_block_ += value;
    if (num_values_ % num_blocks_ == 0) {
      if (block_averages_ == NULL) {
        block_averages_ = std::make_shared<Accumulator>();
      }
      block_averages_->accumulate(sum_block_/static_cast<double>(num_blocks_));
      sum_block_ = 0.;
    }
  }
}

void Accumulator::reset() {
  sum_ = 0;
  num_values_ = 0;
  sum_squared_ = 0;
  set_block();
  sum_block_ = 0;
  max_ = -NEAR_INFINITY;
  min_ = NEAR_INFINITY;
  set_moments();
}

double Accumulator::average() const {
  double av = 0;
  if (num_values_ > 0) {
    av = static_cast<double>(sum_ / static_cast<double>(num_values_));
  }
  return av;
}

double Accumulator::stdev() const {
  double stdev = 0;
  if (num_values_ > 1) {
    const double fluct = sum_squared_/static_cast<double>(num_values_)
                       - pow(average(), 2);
    if (fluct > 0.) {
      stdev = sqrt(fluct*num_values_ / static_cast<double>(num_values_ - 1));
    }
  }
  return stdev;
}

double Accumulator::block_stdev() const {
  if (num_blocks_ != 0) {
    if (block_averages_ != NULL) {
      if (block_averages_->num_values() > 1) {
        return block_averages_->std()/
               sqrt(static_cast<int>(block_averages_->num_values()));
      }
    }
  }
  return 0;
}

double Accumulator::stdev_of_av() const {
  return std()/sqrt(static_cast<double>(num_values()));
}

void Accumulator::set_block(const long long num_block) {
  num_blocks_ = num_block;
}

void Accumulator::set_moments(const int num_moments) {
  val_moment_.resize(num_moments);
}

void Accumulator::serialize(std::ostream& ostr) const {
  feasst_serialize_version(1, ostr);
  feasst_serialize(num_values_, ostr);
  feasst_serialize(sum_, ostr);
  feasst_serialize(sum_squared_, ostr);
  feasst_serialize(max_, ostr);
  feasst_serialize(min_, ostr);
  feasst_serialize(val_moment_, ostr);
  feasst_serialize(num_blocks_, ostr);
  feasst_serialize(sum_block_, ostr);
  feasst_serialize(block_averages_, ostr);
}

Accumulator::Accumulator(std::istream& istr) {
  feasst_deserialize_version(istr);
  feasst_deserialize(&num_values_, istr);
  feasst_deserialize(&sum_, istr);
  feasst_deserialize(&sum_squared_, istr);
  feasst_deserialize(&max_, istr);
  feasst_deserialize(&min_, istr);
  feasst_deserialize(&val_moment_, istr);
  feasst_deserialize(&num_blocks_, istr);
  feasst_deserialize(&sum_block_, istr);
  // feasst_deserialize(block_averages_, istr);
  // HWH for unknown reasons, the above does not work.
  int existing;
  istr >> existing;
  if (existing != 0) {
    block_averages_ = std::make_shared<Accumulator>(istr);
  }
}

}  // namespace feasst

