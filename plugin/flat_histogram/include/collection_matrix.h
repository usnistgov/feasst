
#ifndef FEASST_FLAT_HISTOGRAM_BIAS_COLLECTION_MATRIX_H_
#define FEASST_FLAT_HISTOGRAM_BIAS_COLLECTION_MATRIX_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "math/include/accumulator.h"
#include "flat_histogram/include/ln_probability.h"

namespace feasst {

/**
  The collection matrix is triple banded when the macrostate can only increase
  or decrease by a single "bin" in the macrostate order parameter.
  For example, one cannot use a triple banded collection matrix if insertions of
  both single and pairs of particles are attempted. Single or pair only would be
  fine.

  The collection matrix was slightly reformulated to not require the
  middle band by making the elements intensive.
  Thus, instead of storing C0, C1 and C2, we store P_down=C0/(C0+C1+C2) and
  P_up=C2/(C0+C1+C2).
  Note that C0+C1+C2 is simply the number of trials starting from that
  macrostate.
  Thus, this new, intensive collection matrix is now double banded.

  The first index of the matrix is the macrostate.
  The second is the state change as follows:
    - 0: macrostate decrease
    - 1: macrostate increase
 */
class CollectionMatrix {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - delta_ln_prob_guess: if the CollectionMatrix lacks transitions to compute
      a delta_ln_prob, use this value instead (default: 0).
      This guess can affect the initial convergence.
      When a simulation starts at low macrostate, a negative number encourages
      acceptance of macrostates that have not yet been visited.
      And for high macrostate, a positive number may help.
      If the value is too large, a very low probability trial may be accepted,
      creating a large free energy difference that will make it difficult to
      sample the reverse transition.
      If a simulation appears "stuck" in a subset of the macrostate range,
      then the guess may be too large.
    - visits_per_delta_ln_prob_boost: In the case where a guess is needed,
      if "n" attempts to increase the macrostate, decrease delta_ln_prob_guess
      by 0.01 per "n".
      Alternatively, if "n" attempts to decrease macrostate, increase
      delta_ln_prob_guess by 0.01 per "n".
      If visits_per_delta_ln_prob_boost is -1, do nothing (default: -1).
    - exp_for_boost: ratio of neighboring probabilities must be within
      10^exp or 10^-exp for boosting.
      If -1, ignore (default: 2).
   */
  explicit CollectionMatrix(argtype args = argtype());
  explicit CollectionMatrix(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Construct from a matrix.
  explicit CollectionMatrix(
    const std::vector<std::vector<Accumulator> >& matrix);

  /// Construct from a series of single-state collection matricies.
  explicit CollectionMatrix(
    const std::vector<std::vector<std::vector<Accumulator> > >& data);

  void resize(const int num_macrostates);

  /// Add value for a given macrostate and state change.
  void increment(const int macro, const int state_change, const double add);

  /// Set values
  void set(const int macro, const std::vector<Accumulator>& values);

  /// Update the ln_prob according to the collection matrix.
  void compute_ln_prob(LnProbability * ln_prob,
    /// optionaly compute the ln_prob from a block (if != -1).
    const int block = -1) const;

  /// Return the matrix
  const std::vector<std::vector<Accumulator> >& matrix() const { return matrix_; }

//  /// Return the standard deviation of the change in the ln_prob relative to
//  /// the bin below.
//  double delta_ln_prob_stdev(const int bin, const int block) const;

  std::string write_per_bin(const int bin) const;
  std::string write_per_bin_header() const;

  int min_blocks() const;
  std::vector<LnProbability> ln_prob_blocks() const;

  bool is_equal(const CollectionMatrix& colmat,
      const double tolerance) const;

  void serialize(std::ostream& ostr) const;

  explicit CollectionMatrix(std::istream& istr);

  //@}
 private:
  double delta_ln_prob_guess_;
  int visits_per_delta_ln_prob_boost_;
  double exp_for_boost_;
  std::vector<std::vector<Accumulator> > matrix_;

  int visits_(const int macro, const int block, const bool lower) const;
};

inline std::shared_ptr<CollectionMatrix> MakeCollectionMatrix(
    argtype args = argtype()) {
  return std::make_shared<CollectionMatrix>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_BIAS_COLLECTION_MATRIX_H_
