
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

CollectionMatrix::CollectionMatrix(argtype * args) {
  delta_ln_prob_guess_ = dble("delta_ln_prob_guess", args, 0);
  visits_per_delta_ln_prob_boost_ = integer("visits_per_delta_ln_prob_boost",
    args, -1);
  exp_for_boost_ = dble("exp_for_boost", args, 2);
}
CollectionMatrix::CollectionMatrix(argtype args)
  : CollectionMatrix(&args) {
  FEASST_CHECK_ALL_USED(args);
}

int CollectionMatrix::visits_(const int macro, const int block, const bool lower) const {
  if (block == -1) {
    if (lower) {
      return matrix_[macro][1].num_values();
    } else {
      return matrix_[macro][0].num_values();
    }
  } else {
    if (lower) {
      return matrix_[macro][1].blocks()[0].size();
    } else {
      return matrix_[macro][0].blocks()[0].size();
    }
  }
  FATAL("unrecognized");
}

void CollectionMatrix::compute_ln_prob(
    LnProbability * ln_prob,
    const int block) const {
  ASSERT(ln_prob->size() > 0, "error");
  ln_prob->set_value(0, 0.);
  for (int macro = 1; macro < ln_prob->size(); ++macro) {
    const double ln_prob_previous = ln_prob->value(macro - 1);
//    if (block == -1) INFO("ln_prob_previous " << ln_prob_previous);
    const int vis_up = visits_(macro - 1, block, true);
    const int vis_down = visits_(macro, block, false);
    if (vis_up == 0 || vis_down == 0) {
      double delta_ln_prob = delta_ln_prob_guess_;
      if (visits_per_delta_ln_prob_boost_ > 0) {
        if (vis_up == 0) {
          //delta_ln_prob += 0.01*vis_down/visits_per_delta_ln_prob_boost_;
          bool boost = true;
          if (exp_for_boost_ > 0 && macro < ln_prob->size() - 1) {
            const double ratio = matrix_[macro][0].average()/
                                 matrix_[macro+1][0].average();
            if (ratio > pow(10, exp_for_boost_) || ratio < pow(10, -exp_for_boost_)) {
              boost = false;
              DEBUG("macro " << macro << " ratio " << ratio);
            }
          }
          if (boost) {
            delta_ln_prob -= 0.01*vis_up/visits_per_delta_ln_prob_boost_;
          }
        } else if (vis_down == 0) {
          bool boost = true;
          if (exp_for_boost_ > 0 && macro > 1) {
            const double ratio = matrix_[macro-1][1].average()/
                                 matrix_[macro-2][1].average();
            if (ratio > pow(10, exp_for_boost_) || ratio < pow(10, -exp_for_boost_)) {
              boost = false;
              DEBUG("macro " << macro << " ratio " << ratio);
            }
          }
          if (boost) {
            delta_ln_prob -= 0.01*vis_up/visits_per_delta_ln_prob_boost_;
          }
        }
      }
      ln_prob->set_value(macro, ln_prob_previous + delta_ln_prob);
    } else {
      double prob_decrease;
      if (block == -1) {
        prob_decrease = matrix_[macro][0].average();
      } else {
        prob_decrease = matrix_[macro][0].blocks()[0][block];
      }
//      if (block == -1) INFO("prob_decrease " << prob_decrease);
      if (prob_decrease == 0) {
        ln_prob->set_value(macro, ln_prob_previous);
      } else {
        double prob_increase;
        if (block == -1) {
          prob_increase = matrix_[macro - 1][1].average();
        } else {
          prob_increase = matrix_[macro - 1][1].blocks()[0][block];
        }
//        if (block == -1) INFO("prob_increase " << prob_increase);
        if (prob_increase == 0) {
          ln_prob->set_value(macro, ln_prob_previous);
        } else {
          const double ln_prob_new = ln_prob_previous
                                   + log(prob_increase/prob_decrease);
//          if (block == -1) INFO("ln_prob_new " << ln_prob_new);
          if (std::isnan(ln_prob_new) || std::isinf(ln_prob_new)) {
            ln_prob->set_value(macro, ln_prob_previous);
          } else {
            ln_prob->set_value(macro, ln_prob_new);
          }
        }
      }
    }
  }
  ln_prob->normalize();
}

void CollectionMatrix::increment(
    const int row,
    const int column,
    const double inc) {
  DEBUG("row " << row << " column " << column << " size " << matrix_.size());
  ASSERT(row < static_cast<int>(matrix_.size()),
    "row: " << row << "  size: " << matrix_.size());
  ASSERT(column < static_cast<int>(matrix_[row].size()),
    "row: " << column << "  size: " << matrix_[row].size());
  matrix_[row][column].accumulate(inc);
}

void CollectionMatrix::serialize(std::ostream& ostr) const {
  feasst_serialize_version(2468, ostr);
  feasst_serialize(delta_ln_prob_guess_, ostr);
  feasst_serialize(visits_per_delta_ln_prob_boost_, ostr);
  feasst_serialize(exp_for_boost_, ostr);
  feasst_serialize_fstobj(matrix_, ostr);
}

CollectionMatrix::CollectionMatrix(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2468, "unrecognized verison: " << version);
  feasst_deserialize(&delta_ln_prob_guess_, istr);
  feasst_deserialize(&visits_per_delta_ln_prob_boost_, istr);
  feasst_deserialize(&exp_for_boost_, istr);
  feasst_deserialize_fstobj(&matrix_, istr);
}

bool CollectionMatrix::is_equal(
    const CollectionMatrix& colmat,
    const double tolerance) const {
  for (int row = 0; row < static_cast<int>(matrix_.size()); ++row) {
    for (int col = 0; col < static_cast<int>(matrix_[row].size()); ++col) {
      if (!matrix_[row][col].is_equal(colmat.matrix_[row][col], tolerance)) {
        return false;
      }
    }
  }
  return true;
}

CollectionMatrix::CollectionMatrix(
    const std::vector<std::vector<Accumulator> >& matrix) : CollectionMatrix() {
  matrix_ = matrix;
}
CollectionMatrix::CollectionMatrix(
    const std::vector<std::vector<std::vector<Accumulator> > >& data) :
    CollectionMatrix() {
  for (const std::vector<std::vector<Accumulator> >& dat : data) {
    ASSERT(dat.size() == 1, "only single states accepted for NVT+W");
    matrix_.push_back(dat[0]);
  }
}

std::string CollectionMatrix::write_per_bin_header() const {
  std::stringstream ss;
//  const int block = chosen_block();
//  if (block != -1) {
//    for (int i = 0; i < static_cast<int>(blocks()[block].size()); ++i) {
//      ss << "ln_prob" << i << ",";
//    }
//    ss << "delta_ln_prob_stdev,";
//  }
  //INFO(min_blocks());
  for (int i = 0; i < min_blocks(); ++i) {
    ss << "ln_prob" << i << ",";
  }
  ss << "delta_ln_prob_stdev,P_down,P_up,P_down_block_std,P_up_block_std";
  return ss.str();
}

std::string CollectionMatrix::write_per_bin(const int bin) const {
  std::stringstream ss;
  std::vector<LnProbability> ln_probs = ln_prob_blocks();
  Accumulator delta_ln_p;
  for (const auto& ln_prob : ln_probs) {
    ss << ln_prob.value(bin) << ",";
    delta_ln_p.accumulate(ln_prob.delta(bin));
  }
  ss << delta_ln_p.stdev_of_av() << ",";
  ss << MAX_PRECISION << matrix()[bin][0].average() << ","
     << matrix()[bin][1].average() << ","
     << matrix()[bin][0].block_stdev() << ","
     << matrix()[bin][1].block_stdev();
  return ss.str();
}

void CollectionMatrix::resize(const int num_macrostates) {
  feasst::resize(num_macrostates, 2, &matrix_);
  for (auto& mat1 : matrix_) {
    for (auto& mat2 : mat1) {
      mat2.reset();
    }
  }
}

void CollectionMatrix::set(const int macro, const std::vector<Accumulator>& values) {
  matrix_[macro] = values;
}

int CollectionMatrix::min_blocks() const {
  bool found = false;
  int min = 1e9;
  for (const auto& mat1 : matrix_) {
    for (const auto& mat2 : mat1) {
      if (mat2.block_size().size() > 0) {
        if (mat2.blocks()[0].size() > 0) {
          const int num = static_cast<int>(mat2.blocks()[0].size());
          DEBUG("num " << num);
          if (num < min) {
            min = num;
            found = true;
          }
        }
      }
    }
  }
  if (found) {
    return min;
  } else {
    return 0;
  }
}

std::vector<LnProbability> CollectionMatrix::ln_prob_blocks() const {
  std::vector<LnProbability> ln_probs;
  LnProbability lnpi;
  lnpi.resize(matrix().size());
  for (int block = 0; block < min_blocks(); ++block) {
    compute_ln_prob(&lnpi, block);
    ln_probs.push_back(lnpi);
  }
  return ln_probs;
}

}  // namespace feasst
