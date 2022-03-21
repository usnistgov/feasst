
#include <cmath>
#include <algorithm>
#include "utils/include/utils.h"  // is_equal
#include "utils/include/serialize.h"
#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "math/include/accumulator.h"
#include "flat_histogram/include/collection_matrix.h"

namespace feasst {

TripleBandedCollectionMatrix::TripleBandedCollectionMatrix(argtype * args) {
  max_block_operations_ = integer("max_block_operations", args, 5);
  updates_ = 0;
}
TripleBandedCollectionMatrix::TripleBandedCollectionMatrix(argtype args)
  : TripleBandedCollectionMatrix(&args) {
  check_all_used(args);
}

void TripleBandedCollectionMatrix::compute_ln_prob(
    LnProbability * ln_prob) const {
  ln_prob->set_value(0, 0.);
  for (int macro = 1; macro < ln_prob->size(); ++macro) {
    const double ln_prob_previous = ln_prob->value(macro - 1);
    TRACE("ln_prob_previous " << ln_prob_previous);
    const double collect_sum_previous = sum_(macro - 1);
    if (collect_sum_previous == 0) {
      ln_prob->set_value(macro, ln_prob_previous);
    } else {
      const double collect_sum = sum_(macro);
      if (collect_sum == 0) {
        ln_prob->set_value(macro, ln_prob_previous);
      } else {
        const double prob_decrease = matrix_[macro][0]/collect_sum;
        TRACE("prob_decrease " << prob_decrease);
        if (prob_decrease == 0) {
          ln_prob->set_value(macro, ln_prob_previous);
        } else {
          const double prob_increase = matrix_[macro - 1][2]/collect_sum_previous;
          TRACE("prob_increase " << prob_increase);
          if (prob_increase == 0) {
            ln_prob->set_value(macro, ln_prob_previous);
          } else {
            const double ln_prob_new = ln_prob_previous
                                     + log(prob_increase/prob_decrease);
            TRACE("ln_prob_new " << ln_prob_new);
            ASSERT(!std::isnan(ln_prob_new), "error");
            ASSERT(std::abs(ln_prob_new) < NEAR_INFINITY, "ln_prob_new:" << ln_prob_new);
            ln_prob->set_value(macro, ln_prob_new);
          }
        }
      }
    }
  }
  ln_prob->normalize();
}

void TripleBandedCollectionMatrix::init_cur_(const int exp) {
  cur_block_[exp] = TripleBandedCollectionMatrix();
  cur_block_[exp].block_ = true;
  cur_block_[exp].resize(matrix().size());
}

void TripleBandedCollectionMatrix::increment(
    const int row,
    const int column,
    const double inc) {
  DEBUG("row " << row << " column " << column << " size " << matrix_.size());
  DEBUG(matrix_[row][column] << " + " << inc);
  ASSERT(row < static_cast<int>(matrix_.size()),
    "row: " << row << "  size: " << matrix_.size());
  ASSERT(column < static_cast<int>(matrix_[row].size()),
    "row: " << column << "  size: " << matrix_[row].size());
  matrix_[row][column] += inc;
  DEBUG(feasst_str(matrix_[row]));

  if (block_) return;

  if (max_block_operations_ > 0) {
    ++updates_;
    // initialize
    if (blocks_.size() == 0) {
      blocks_.resize(max_block_operations_);
      cur_block_.resize(max_block_operations_);
      max_block_updates_.resize(max_block_operations_);
      block_updates_.resize(max_block_operations_, 0);
      for (int exp = 0; exp < max_block_operations_; ++exp) {
        init_cur_(exp);
        max_block_updates_[exp] = std::pow(2, exp);
      }
    }

    // check for new size
    const long double new_block_size = 2*max_block_updates_.back();
    if (std::abs(std::fmod(updates_, new_block_size)) < 0.1) {
      max_block_updates_.erase(max_block_updates_.begin());
      max_block_updates_.push_back(new_block_size);
      block_updates_.erase(block_updates_.begin());
      block_updates_.push_back(updates_ - 1);
      blocks_.erase(blocks_.begin());
      blocks_.push_back(std::vector<TripleBandedCollectionMatrix>());
      cur_block_.erase(cur_block_.begin());
      //cur_block_.push_back(*this);
      cur_block_.push_back(TripleBandedCollectionMatrix());
      init_cur_(max_block_operations_ - 1);
      cur_block_.back().matrix_ = matrix_;
    }

    for (int exp = 0; exp < max_block_operations_; ++exp) {
      cur_block_[exp].increment(row, column, inc);
      ++block_updates_[exp];
      if (block_updates_[exp] == max_block_updates_[exp]) {
        block_updates_[exp] = 0;
        blocks_[exp].push_back(cur_block_[exp]);
        init_cur_(exp);
      }
    }
  }
}

double TripleBandedCollectionMatrix::sum_(const int macro) const {
  return std::accumulate(matrix_[macro].begin(), matrix_[macro].end(), 0.);
}

void TripleBandedCollectionMatrix::serialize(std::ostream& ostr) const {
  feasst_serialize_version(2468, ostr);
  feasst_serialize(matrix_, ostr);
  feasst_serialize(block_, ostr);
  feasst_serialize(max_block_operations_, ostr);
  feasst_serialize(updates_, ostr);
  feasst_serialize(block_updates_, ostr);
  feasst_serialize(max_block_updates_, ostr);
  feasst_serialize_fstobj(blocks_, ostr);
  feasst_serialize_fstobj(cur_block_, ostr);
}

TripleBandedCollectionMatrix::TripleBandedCollectionMatrix(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2468, "unrecognized verison: " << version);
  feasst_deserialize(&matrix_, istr);
  feasst_deserialize(&block_, istr);
  feasst_deserialize(&max_block_operations_, istr);
  feasst_deserialize(&updates_, istr);
  feasst_deserialize(&block_updates_, istr);
  feasst_deserialize(&max_block_updates_, istr);
  feasst_deserialize_fstobj(&blocks_, istr);
  feasst_deserialize_fstobj(&cur_block_, istr);
}

bool TripleBandedCollectionMatrix::is_equal(
    const TripleBandedCollectionMatrix& colmat,
    const double tolerance) const {
  if (!feasst::is_equal(matrix_, colmat.matrix_, tolerance)) {
    INFO("colmat not equal " << feasst_str(matrix_, true));
    INFO("colmat not equal " << feasst_str(colmat.matrix_, true));
    return false;
  }
  return true;
}

TripleBandedCollectionMatrix::TripleBandedCollectionMatrix(
    const std::vector<std::vector<double> >& matrix) {
  matrix_ = matrix;
}
TripleBandedCollectionMatrix::TripleBandedCollectionMatrix(
    const std::vector<std::vector<std::vector<double> > >& data) {
  for (const std::vector<std::vector<double> >& dat : data) {
    ASSERT(dat.size() == 1, "only single states accepted for NVT+W");
    matrix_.push_back(dat[0]);
  }
}

int TripleBandedCollectionMatrix::chosen_block_() const {
  int bl = static_cast<int>(blocks().size()) - 1;
  while (bl >= 0) {
    if (static_cast<int>(blocks()[bl].size()) >= 10) {
      return bl;
    }
    --bl;
  }
  return bl;
}

std::string TripleBandedCollectionMatrix::write_per_bin_header() const {
  std::stringstream ss;
  const int chosen_block = chosen_block_();
  if (chosen_block != -1) {
    for (int i = 0; i < static_cast<int>(blocks()[chosen_block].size()); ++i) {
      ss << "ln_prob" << i << ",";
    }
    ss << "delta_ln_prob_stdev,";
  }
  ss << "c0,c1,c2,";
  return ss.str();
}

std::string TripleBandedCollectionMatrix::write_per_bin(const int bin) const {
  std::stringstream ss;
  // compute and print the macrostates of each block
  // and calculate the standard deviation of the average of the blocks
  const int chosen_block = chosen_block_();
  if (chosen_block != -1) {
    Accumulator acc_block;
    for (const TripleBandedCollectionMatrix& cm : blocks()[chosen_block]) {
      LnProbability lnp;
      lnp.resize(matrix().size());
      cm.compute_ln_prob(&lnp);
      ss << std::setprecision(5) << lnp.value(bin) << ",";
      if (bin == 0) {
        acc_block.accumulate(0);
      } else {
        acc_block.accumulate(lnp.value(bin) - lnp.value(bin-1));
      }
    }
    ss << std::setprecision(4) << acc_block.stdev_of_av() << ",";
  }
  ss << MAX_PRECISION << matrix()[bin][0] << ","
     << matrix()[bin][1] << ","
     << matrix()[bin][2] << ",";
//    ss << " " << cols[0] << " "
//      << cols[1] << " "
//      << cols[2];
  return ss.str();
}

void TripleBandedCollectionMatrix::resize(const int num_macrostates) {
  feasst::resize(num_macrostates, 3, &matrix_);
  feasst::fill(0., &matrix_);
}

void TripleBandedCollectionMatrix::set(const int macro, const std::vector<double>& values) {
  matrix_[macro] = values;
}

}  // namespace feasst
