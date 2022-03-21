
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

  Blocks and error bars are obtained via the blocking method.
  Blocks are created on the fly by sizes of base two.
  The default blocks and stdev are from the largest block sizes with more than
  ten blocks.
  Error bars on grand canonical ensemble averages may then by obtained by using
  the macrostate distribution from each of these blocks.
 */
class TripleBandedCollectionMatrix {
 public:
  /**
    args:
    - max_block_operations: maximum number of blocking operations (default: 5).
   */
  explicit TripleBandedCollectionMatrix(argtype args = argtype());
  explicit TripleBandedCollectionMatrix(argtype * args);

  /// Construct from a matrix.
  explicit TripleBandedCollectionMatrix(
    const std::vector<std::vector<double> >& matrix);

  /// Construct from a series of single-state collection matricies.
  explicit TripleBandedCollectionMatrix(
    const std::vector<std::vector<std::vector<double> > >& data);

  void resize(const int num_macrostates);

  /// Add value for a given macrostate and state change.
  void increment(const int macro, const int state_change, const double add);

  /// Set values
  void set(const int macro, const std::vector<double>& values);

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
  int max_block_operations_;
  long double updates_;
  std::vector<double> block_updates_;
  std::vector<double> max_block_updates_;
  std::vector<std::vector<TripleBandedCollectionMatrix> > blocks_;
  std::vector<TripleBandedCollectionMatrix> cur_block_;

  void init_cur_(const int exp);
  double sum_(const int macro) const;

  // find the largest block size with >= 10 blocks.
  // if none, return -1
  int chosen_block_() const;
};

inline std::shared_ptr<TripleBandedCollectionMatrix> MakeTripleBandedCollectionMatrix(
    argtype args = argtype()) {
  return std::make_shared<TripleBandedCollectionMatrix>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_BIAS_COLLECTION_MATRIX_H_
