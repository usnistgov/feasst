#include <cmath>
#include <numeric>
#include <sstream>
#include "utils/include/io.h"
#include "utils/include/serialize.h"
#include "math/include/accumulator.h"
#include "math/include/constants.h"

namespace feasst {

Accumulator::Accumulator(argtype args) {
  reset();
  set_block_size(integer("block_size", &args, 1e5));
  check_all_used(args);
}

// Accumulator::Accumulator(const long long num_values, const long double sum,
//   const long double sum_squared) {
//   reset();
//   num_values_ = num_values;
//   sum_ = sum;
//   sum_squared_ = sum_squared;
// }

void Accumulator::accumulate(double value) {
  last_value_ = value;
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
  if (block_size_ != 0) {
    sum_block_ += value;
    if (std::abs(std::fmod(num_values_, block_size_)) < 0.1) {
      if (block_averages_ == NULL) {
        block_averages_ = std::make_shared<Accumulator>();
      }
      block_averages_->accumulate(sum_block_/block_size_);
      sum_block_ = 0.;
    }
  }
}

void Accumulator::reset() {
  sum_ = 0;
  num_values_ = 0;
  sum_squared_ = 0;
  set_block_size();
  sum_block_ = 0;
  max_ = -NEAR_INFINITY;
  min_ = NEAR_INFINITY;
  set_moments();
}

double Accumulator::average() const {
  double av = 0;
  if (num_values_ > 0) {
    av = static_cast<double>(sum_/num_values_);
  }
  return av;
}

double Accumulator::stdev() const {
  double stdev = 0;
  if (num_values_ > 1) {
    const double fluct = sum_squared_/num_values_
                       - std::pow(average(), 2);
    if (fluct > 0.) {
      stdev = sqrt(fluct*num_values_/(num_values_ - 1));
    }
  }
  return stdev;
}

double Accumulator::block_stdev() const {
  if (block_size_ != 0) {
    if (block_averages_ != NULL) {
      if (block_averages_->num_values() > 1) {
        return block_averages_->std()/
               std::sqrt(block_averages_->num_values());
      }
    }
  }
  return 0;
}

double Accumulator::stdev_of_av() const {
  return std()/std::sqrt(num_values());
}

void Accumulator::set_moments(const int num_moments) {
  val_moment_.resize(num_moments);
}

std::string Accumulator::status_header() const {
  std::stringstream ss;
  ss << "average,stdev,block_stdev,n,";
  for (int i = 0; i < static_cast<int>(val_moment_.size()); ++i) {
    ss << "moment" << i << ",";
  }
  return ss.str();
}

std::string Accumulator::status() const {
  std::stringstream ss;
  ss << MAX_PRECISION << average() << ",";
  ss << stdev() << ",";
  ss << block_stdev() << ",";
  ss << num_values_ << ",";
  for (const long double moment : val_moment_) {
    ss << moment << ",";
  }
  return ss.str();
}

std::string Accumulator::str() const {
  std::stringstream ss;
  ss << status_header() << std::endl << status();
  return ss.str();
}

void Accumulator::serialize(std::ostream& ostr) const {
  feasst_serialize_version(8773, ostr);
  feasst_serialize(num_values_, ostr);
  feasst_serialize(sum_, ostr);
  feasst_serialize(sum_squared_, ostr);
  feasst_serialize(max_, ostr);
  feasst_serialize(min_, ostr);
  feasst_serialize(val_moment_, ostr);
  feasst_serialize(block_size_, ostr);
  feasst_serialize(sum_block_, ostr);
  feasst_serialize(block_averages_, ostr);
}

Accumulator::Accumulator(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8773, "unrecognized verison: " << version);
  feasst_deserialize(&num_values_, istr);
  feasst_deserialize(&sum_, istr);
  feasst_deserialize(&sum_squared_, istr);
  feasst_deserialize(&max_, istr);
  feasst_deserialize(&min_, istr);
  feasst_deserialize(&val_moment_, istr);
  feasst_deserialize(&block_size_, istr);
  feasst_deserialize(&sum_block_, istr);
  // feasst_deserialize(block_averages_, istr);
  // HWH for unknown reasons, the above does not work.
  int existing;
  istr >> existing;
  if (existing != 0) {
    block_averages_ = std::make_shared<Accumulator>(istr);
  }
}

double Accumulator::last_value() const {
  ASSERT(num_values_ > 0, "no values accumulated");
  return last_value_;
}

bool Accumulator::is_equivalent(const Accumulator& accum,
                                const double t_factor,
                                const bool verbose) const {
  const double diff = accum.average() - average();
  const double stdev = std::sqrt(
    pow(accum.block_stdev(), 2)/accum.block_averages_->num_values() +
    pow(block_stdev(), 2)/block_averages_->num_values());
  if (std::abs(diff) > t_factor*stdev) {
    if (verbose) {
      std::cout << str() << " " << accum.str() << std::endl;
      std::cout << "diff: " << diff << std::endl;
      std::cout << "stdev: " << stdev << std::endl;
      std::cout << "t_factor: " << t_factor << std::endl;
    }
    return false;
  }
  return true;
}

}  // namespace feasst
