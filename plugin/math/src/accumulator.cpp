#include <cmath>
#include <numeric>
#include <sstream>
#include "utils/include/utils.h"
#include "utils/include/io.h"
#include "utils/include/serialize.h"
#include "math/include/accumulator.h"
#include "math/include/constants.h"

namespace feasst {

Accumulator::Accumulator(argtype * args) {
  max_block_operations_ = integer("max_block_operations", args, 6);
  set_moments_(integer("num_moments", args, 5));
  reset();
}

Accumulator::Accumulator(argtype args) : Accumulator(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void Accumulator::accumulate(double value) {
  last_value_ = value;
  long double valmo = 1.;
  for (int mo = 0; mo < static_cast<int>(val_moment_.size()); ++mo) {
    val_moment_[mo] += valmo;
    valmo *= value;
  }
  if (max_ < value) max_ = value;
  if (min_ > value) min_ = value;

  // check if its time to make a new block size
  if (max_block_operations_ > 0) {
    const double new_block_size = block_power_*block_size_.back();
    //INFO("new_block_size " << new_block_size);
    if (std::abs(std::fmod(num_values(), new_block_size)) < 0.1) {
      block_size_.erase(block_size_.begin());
      block_size_.push_back(new_block_size);
      sum_block_.erase(sum_block_.begin());
      sum_block_.push_back(sum());
      block_averages_.erase(block_averages_.begin());
      block_averages_.push_back(MakeAccumulator({{"max_block_operations", "0"},
        {"num_moments", feasst::str(num_moments())}}));
      blocks_.erase(blocks_.begin());
      blocks_.push_back(std::vector<double>());
    }

    // accumulate block averages
    for (int bop = 0; bop < max_block_operations_; ++bop) {
      //INFO("bop " << bop << " size " << block_size_[bop]);
      sum_block_[bop] += value;
      if (std::abs(std::fmod(num_values(), block_size_[bop])) < 0.1) {
        ASSERT(block_averages_[bop], "er");
        const double new_av = sum_block_[bop]/block_size_[bop];
        block_averages_[bop]->accumulate(new_av);
        blocks_[bop].push_back(new_av);
        sum_block_[bop] = 0.0L;
      }
    }
  }
}

void Accumulator::reset() {
  max_ = -NEAR_INFINITY;
  min_ = NEAR_INFINITY;
  for (long double& mom : val_moment_) {
    mom = 0.0L;
  }
  if (static_cast<int>(block_averages_.size()) > 0) {
    block_averages_.clear();
    block_size_.clear();
    blocks_.clear();
  }
  block_averages_.resize(max_block_operations_);
  block_size_.resize(max_block_operations_);
  blocks_.resize(max_block_operations_);
  sum_block_.resize(max_block_operations_, 0.);
  for (int bop = 0; bop < max_block_operations_; ++bop) {
    block_size_[bop] = std::pow(block_power_, bop);
    block_averages_[bop] = MakeAccumulator({{"max_block_operations", "0"},
      {"num_moments", feasst::str(num_moments())}});
    blocks_[bop].clear();
  }
}

double Accumulator::average() const {
  double av = 0;
  if (num_values() > 0) {
    av = static_cast<double>(sum()/num_values());
  }
  return av;
}

double Accumulator::stdev() const {
  double stdev = 0;
  if (num_values() > 1) {
    const double fluct = sum_of_squared()/num_values()
                       - std::pow(average(), 2);
    if (fluct > 0.) {
      stdev = sqrt(fluct*num_values()/(num_values() - 1));
    }
  }
  return stdev;
}

double Accumulator::block_stdev(const int num_op, const int min_blocks) const {
  double std = 0.;
  if (num_op == -1) {
    for (std::shared_ptr<Accumulator> acc : block_averages()) {
      if (acc->num_values() > min_blocks) {
        const double acc_std = acc->stdev_of_av();
        if (acc_std > std) {
          std = acc_std;
        }
      }
    }
  } else {
    if (block_averages_[num_op] != NULL) {
      if (block_averages_[num_op]->num_values() > 1) {
        std =  block_averages_[num_op]->std()/
               std::sqrt(block_averages_[num_op]->num_values());
      }
    }
  }
  return std;
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
  ss << "average,stdev,block_stdev,";
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
  feasst_serialize(max_, ostr);
  feasst_serialize(min_, ostr);
  feasst_serialize(val_moment_, ostr);
  feasst_serialize(max_block_operations_, ostr);
  feasst_serialize(block_size_, ostr);
  feasst_serialize(sum_block_, ostr);
  feasst_serialize(blocks_, ostr);
  feasst_serialize(block_averages_, ostr);
}

Accumulator::Accumulator(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8773, "unrecognized verison: " << version);
  feasst_deserialize(&max_, istr);
  feasst_deserialize(&min_, istr);
  feasst_deserialize(&val_moment_, istr);
  feasst_deserialize(&max_block_operations_, istr);
  feasst_deserialize(&block_size_, istr);
  feasst_deserialize(&sum_block_, istr);
  feasst_deserialize(&blocks_, istr);
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
  ASSERT(num_values() > 0, "no values accumulated");
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
  ASSERT(num_moments() > n, "number of moments: " << num_moments()
    << " must be > n: " << n);
  const double v = num_values();
  if (n == 2) {
    return moment(2)/v - std::pow(moment(1)/v, 2);
    //return std();
  } else if (n == 3) {
    return moment(3)/v - 3.*moment(1)*moment(2)/v/v + 2.*std::pow(moment(1)/v, 3);
  } else if (n == 4) {
    return moment(4)/v-4*moment(1)*moment(3)/v/v+6*moment(1)*moment(1)*moment(2)/v/v/v-3*std::pow(moment(1)/v, 4);
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

bool Accumulator::is_equal(const Accumulator& acc, const double tolerance) const {
  if (max_ != acc.max_) return false;
  if (min_ != acc.min_) return false;
  if (!feasst::is_equal(val_moment_, acc.val_moment_, tolerance)) return false;
  if (block_power_ != acc.block_power_) return false;
  if (!feasst::is_equal(block_size_, acc.block_size_, tolerance)) return false;
  for (int av = 0; av < static_cast<int>(block_averages_.size()); ++av) {
    if (!block_averages_[av]->is_equal(*acc.block_averages_[av], tolerance)) {
      return false;
    }
  }
  if (!feasst::is_equal(blocks_, acc.blocks_, tolerance)) return false;
  return true;
}

}  // namespace feasst
