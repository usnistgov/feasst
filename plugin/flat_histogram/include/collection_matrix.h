
#ifndef FEASST_FLAT_HISTOGRAM_BIAS_COLLECTION_MATRIX_H_
#define FEASST_FLAT_HISTOGRAM_BIAS_COLLECTION_MATRIX_H_

#include <vector>
#include <memory>
#include "flat_histogram/include/ln_probability.h"

namespace feasst {

/**
  The collection matrix is triple banded when the macrostate can only increase
  or decrease by a single "bin" in the macrostate order parameter.

  For example, one cannot use a triple banded collection matrix if insertions of
  both single and pairs of particles are attempted. Single or pair only would be
  fine.

  The first index of the matrix is the macrostate.
  The second is the state change as follows:
    - 0: macrostate decrease
    - 1: macrostate constant
    - 2: macrostate increase

  Beware that the computed probability distribution from this collection matrix
  may contain spurious values at the two most extreme ends.
  These two extreme ends are the two pairs of the two lowest and highest
  indices, respectively.
  For example, at the highest macrostate, a transition to the next highest
  macrostate is immediately rejected and thus the macrostate increase part of
  the collection matrix is zero due to this constraint.
  However, this would not be the case if the higher macrostate values were
  allowed.
  The constraint thus affects the probability at the most extreme values, and
  these values should be discarded.
  There are exceptions such as the case with zero particles.
  In this case, there cannot be less than zero particles so the constraint is
  physically meaningful and the probability at this macrostate is correct.
 */
class TripleBandedCollectionMatrix {
 public:
  TripleBandedCollectionMatrix() {}

  void resize(const int num_macrostates) {
    matrix_.resize(num_macrostates, std::vector<double>(3, 0.)); }

  /// Add value for a given macrostate and state change.
  void increment(const int macro, const int state_change, const double add);

  /// Update the ln_prob according to the collection matrix.
  void compute_ln_prob(LnProbability * ln_prob);

  /// Return the matrix
  std::vector<std::vector<double> > matrix() const { return matrix_; }

  bool is_equal(const TripleBandedCollectionMatrix& colmat,
      const double tolerance) const;

  void serialize(std::ostream& ostr) const;

  TripleBandedCollectionMatrix(std::istream& istr);

 private:
  std::vector<std::vector<double> > matrix_;

  double sum_(const int macro);
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_BIAS_COLLECTION_MATRIX_H_
