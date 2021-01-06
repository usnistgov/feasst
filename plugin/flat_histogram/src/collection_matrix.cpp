
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

void TripleBandedCollectionMatrix::compute_ln_prob(
    LnProbability * ln_prob) {
  ln_prob->set_value(0, 0.);
  for (int macro = 1; macro < ln_prob->size(); ++macro) {
    const double ln_prob_previous = ln_prob->value(macro - 1);
    INFO("ln_prob_previous " << ln_prob_previous);
    const double collect_sum_previous = sum_(macro - 1);
    if (collect_sum_previous == 0) {
      ln_prob->set_value(macro, ln_prob_previous);
    } else {
      const double collect_sum = sum_(macro);
      if (collect_sum == 0) {
        ln_prob->set_value(macro, ln_prob_previous);
      } else {
        const double prob_decrease = matrix_[macro][0]/collect_sum;
        INFO("prob_decrease " << prob_decrease);
        if (prob_decrease == 0) {
          ln_prob->set_value(macro, ln_prob_previous);
        } else {
          const double prob_increase = matrix_[macro - 1][2]/collect_sum_previous;
          INFO("prob_increase " << prob_increase);
          if (prob_increase == 0) {
            ln_prob->set_value(macro, ln_prob_previous);
          } else {
            const double ln_prob_new = ln_prob_previous
                                     + log(prob_increase/prob_decrease);
            INFO("ln_prob_new " << ln_prob_new);
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

void TripleBandedCollectionMatrix::increment(
    const int row,
    const int column,
    const double inc) {
  DEBUG("row " << row << " column " << column << " size " << matrix_.size());
  DEBUG(matrix_[row][column] << " + " << inc);
  matrix_[row][column] += inc;
  DEBUG(feasst_str(matrix_[row]));
}

double TripleBandedCollectionMatrix::sum_(const int macro) {
  return std::accumulate(matrix_[macro].begin(), matrix_[macro].end(), 0.);
}

void TripleBandedCollectionMatrix::serialize(std::ostream& ostr) const {
  feasst_serialize_version(2468, ostr);
  feasst_serialize(matrix_, ostr);
}

TripleBandedCollectionMatrix::TripleBandedCollectionMatrix(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2468, "unrecognized verison: " << version);
  feasst_deserialize(&matrix_, istr);
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

}  // namespace feasst
