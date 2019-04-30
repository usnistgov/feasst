
#include <cmath>
#include <algorithm>
#include "flat_histogram/include/collection_matrix.h"
#include "core/include/utils_io.h"
#include "core/include/utils_math.h"
#include "core/include/debug.h"
#include "core/include/accumulator.h"

namespace feasst {

void TripleBandedCollectionMatrix::compute_ln_prob(
    LnProbabilityDistribution * ln_prob) {
  ln_prob->set_value(0, 0.);
  for (int macro = 1; macro < ln_prob->size(); ++macro) {
    const double ln_prob_previous = ln_prob->value(macro - 1);
    DEBUG("ln_prob_previous " << ln_prob_previous);
    const double collect_sum_previous = sum_(macro - 1);
    if (collect_sum_previous == 0) {
      ln_prob->set_value(macro, ln_prob_previous);
    } else {
      const double collect_sum = sum_(macro);
      if (collect_sum == 0) {
        ln_prob->set_value(macro, ln_prob_previous);
      } else {
        const double prob_decrease = matrix_[macro][0]/collect_sum;
        DEBUG("prob_decrease " << prob_decrease);
        if (prob_decrease == 0) {
          ln_prob->set_value(macro, ln_prob_previous);
        } else {
          const double prob_increase = matrix_[macro - 1][2]/collect_sum_previous;
          DEBUG("prob_increase " << prob_increase);
          if (prob_increase == 0) {
            ln_prob->set_value(macro, ln_prob_previous);
          } else {
            const double ln_prob_new = ln_prob_previous
                                     + log(prob_increase/prob_decrease);
            DEBUG("ln_prob_new " << ln_prob_new);
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
  matrix_[row][column] += inc;
}

double TripleBandedCollectionMatrix::sum_(const int macro) {
  return std::accumulate(matrix_[macro].begin(), matrix_[macro].end(), 0.);
}

void TripleBandedCollectionMatrix::serialize(std::ostream& ostr) const {
  feasst_serialize_version(1, ostr);
  feasst_serialize(matrix_, ostr);
}

TripleBandedCollectionMatrix::TripleBandedCollectionMatrix(std::istream& istr) {
  feasst_deserialize_version(istr);
  feasst_deserialize(&matrix_, istr);
}

}  // namespace feasst
