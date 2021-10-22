#include <cmath>
#include <numeric>
#include <sstream>
#include "utils/include/io.h"
#include "utils/include/serialize.h"
#include "math/include/accumulator.h"
#include "math/include/constants.h"

namespace feasst {

Accumulator::Accumulator(argtype * args) {
  reset();
  num_block_operations_ = integer("num_block_operations", args, 20);
  min_block_size_ = integer("min_block_size", args, 1e5);
  sum_block_.resize(num_block_operations_, 0.);
  block_averages_.resize(num_block_operations_);
  set_moments_(integer("num_moments", args, 4));
}

Accumulator::Accumulator(argtype args) : Accumulator(&args) {
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
  for (int bop = 0; bop < num_block_operations_; ++bop) {
    sum_block_[bop] += value;
    const int block_size = min_block_size_ * std::pow(2, bop);
    if (std::abs(std::fmod(num_values_, block_size)) < 0.1) {
      if (block_averages_[bop] == NULL) {
        block_averages_[bop] = MakeAccumulator({{"num_block_operations", "0"},
          {"num_moments", feasst::str(num_moments())}});
      }
      block_averages_[bop]->accumulate(sum_block_[bop]/block_size);
      sum_block_[bop] = 0.;
    }
  }
}

void Accumulator::reset() {
  num_values_ = 0;
  sum_ = 0;
  sum_squared_ = 0;
  max_ = -NEAR_INFINITY;
  min_ = NEAR_INFINITY;
  for (long double& sum : sum_block_) {
    sum = 0.;
  }
  for (long double& mom : val_moment_) {
    mom = 0.;
  }
  if (static_cast<int>(block_averages_.size()) > 0) {
    block_averages_.clear();
    block_averages_.resize(num_block_operations_);
  }
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

double Accumulator::block_stdev(const int num_op) const {
  if (block_averages_[num_op] != NULL) {
    if (block_averages_[num_op]->num_values() > 1) {
      return block_averages_[num_op]->std()/
             std::sqrt(block_averages_[num_op]->num_values());
    }
  }
  return 0;
}

double Accumulator::stdev_of_av() const {
  return std()/std::sqrt(num_values());
}

void Accumulator::set_moments_(const int num_moments) {
  ASSERT(num_moments > 0, "num_moments: " << num_moments << " >0");
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
  feasst_serialize(num_block_operations_, ostr);
  feasst_serialize(min_block_size_, ostr);
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
  feasst_deserialize(&num_block_operations_, istr);
  feasst_deserialize(&min_block_size_, istr);
  feasst_deserialize(&sum_block_, istr);
  // feasst_deserialize(block_averages_, istr);
  // HWH for unknown reasons, the above does not work.
  int dim1;
  istr >> dim1;
  block_averages_.resize(dim1);
  for (int index = 0; index < dim1; ++index) {
    int existing;
    istr >> existing;
    if (existing != 0) {
      block_averages_[index] = std::make_shared<Accumulator>(istr);
    }
  }
}

double Accumulator::last_value() const {
  ASSERT(num_values_ > 0, "no values accumulated");
  return last_value_;
}

bool Accumulator::is_equivalent(const Accumulator& accum,
                                const double t_factor,
                                const int num_op,
                                const bool verbose) const {
  const double diff = accum.average() - average();
  const double stdev = std::sqrt(
    pow(accum.block_stdev(num_op), 2)/accum.block_averages_[num_op]->num_values() +
    pow(block_stdev(num_op), 2)/block_averages_[num_op]->num_values());
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

double Accumulator::central_moment(const int n) {
  //INFO(val_moment_.size());
  ASSERT(num_moments() >= n, "number of moments: " << num_moments()
    << " must be >= n: " << n);
  const double v = num_values_;
  if (n == 2) {
    return moment(1)/v - std::pow(moment(0)/v, 2);
    //return std();
  } else if (n == 3) {
    return moment(2)/v - 3.*moment(0)*moment(1)/v/v + 2.*std::pow(moment(0)/v, 3);
  } else if (n == 4) {
    return moment(3)/v-4*moment(0)*moment(2)/v/v+6*moment(0)*moment(0)*moment(1)/v/v/v-3*std::pow(moment(0)/v, 4);
  } else {
    FATAL("moment: " << n << " not recognized");
  }
}

double Accumulator::block_std_of_std(const int op) const {
  if (block_averages_[op]) {
    const double num_val = block_averages_[op]->num_values();
    if (num_val > 2) {
      return std::pow(2*block_averages_[op]->central_moment(4)/
                      std::pow(num_val - 1, 3), 1./4.);
    }
  }
  return 0.;
}

}  // namespace feasst
