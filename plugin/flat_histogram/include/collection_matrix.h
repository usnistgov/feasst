
#ifndef FEASST_FLAT_HISTOGRAM_BIAS_COLLECTION_MATRIX_H_
#define FEASST_FLAT_HISTOGRAM_BIAS_COLLECTION_MATRIX_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
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
  This is especially true in cases where windows with different moves sets are
  spliced.
 */
class TripleBandedCollectionMatrix {
 public:
  /**
    args:
    - min_block_size: minimum number of increments for a block (default: 1e6).
    - num_block_operations: maximum exponent that determines maxum number of increments for
      a block (default: 20).
   */
  explicit TripleBandedCollectionMatrix(argtype args = argtype());
  explicit TripleBandedCollectionMatrix(argtype * args);

  /// Construct from a series of single-state collection matricies (NVT+W).
  explicit TripleBandedCollectionMatrix(
    const std::vector<std::vector<std::vector<double> > >& data);

  void resize(const int num_macrostates) {
    matrix_.resize(num_macrostates, std::vector<double>(3, 0.)); }

  /// Add value for a given macrostate and state change.
  void increment(const int macro, const int state_change, const double add);

  /// Update the ln_prob according to the collection matrix.
  void compute_ln_prob(LnProbability * ln_prob) const;

  /// Return the matrix
  const std::vector<std::vector<double> >& matrix() const { return matrix_; }

  std::string write_per_bin(const int bin) const;
  std::string write_per_bin_header() const;

  bool is_equal(const TripleBandedCollectionMatrix& colmat,
      const double tolerance) const;

  void serialize(std::ostream& ostr) const;

  const std::vector<std::vector<TripleBandedCollectionMatrix> >& blocks() const {
    return blocks_; }

  explicit TripleBandedCollectionMatrix(std::istream& istr);

 private:
  std::vector<std::vector<double> > matrix_;

  // blocks
  bool block_ = false;
  int min_block_size_;
  int num_block_operations_;
  std::vector<int> block_updates_;
  std::vector<int> max_block_updates_;
  std::vector<std::vector<TripleBandedCollectionMatrix> > blocks_;
  std::vector<TripleBandedCollectionMatrix> cur_block_;

  void init_cur_(const int exp);

  double sum_(const int macro) const;
};

inline std::shared_ptr<TripleBandedCollectionMatrix> MakeTripleBandedCollectionMatrix(
    argtype args = argtype()) {
  return std::make_shared<TripleBandedCollectionMatrix>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_BIAS_COLLECTION_MATRIX_H_
