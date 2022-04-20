
#ifndef FEASST_FLAT_HISTOGRAM_TRANSITION_MATRIX_H_
#define FEASST_FLAT_HISTOGRAM_TRANSITION_MATRIX_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "flat_histogram/include/bias.h"
#include "flat_histogram/include/collection_matrix.h"

namespace feasst {

/**
  Transition matrix flat histogram bias.

  Direct calculation of liquid–vapor phase equilibria from transition matrix
  Monte Carlo simulation
  https://doi.org/10.1063/1.1572463

  Elucidating the effects of adsorbent flexibility on fluid adsorption using
  simple models and flat-histogram sampling methods.
  https://doi.org/10.1063/1.4884124
 */
class TransitionMatrix : public Bias {
 public:
  /**
    args:
    - min_visits: A sweep is performed when all macrostates are visited by
      another macrostate this number of times (default: 100).
    - average_visits: A sweep is performed when macrostates are visited by
      another macrostate more than this average number of times (default: 0).
    - min_sweeps: Number of sweeps required for completion.
    - reset_sweeps: The 'phase' counter increments from 0 to 1 when
      reset_sweeps are completed (default: -1 [counter will never increment])
    - new_sweep: if set to 1, use new sweep definition of "the minimum number
      of transitions out of each macrostate" (default: 0).
   */
  explicit TransitionMatrix(argtype args = argtype());
  explicit TransitionMatrix(argtype * args);
  void update_or_revert(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted,
    const bool is_endpoint,
    const bool revert) override;

  /// Return the minimum sweeps required for completion.
  int min_sweeps() const { return min_sweeps_; }

  /// Return the reset_sweeps parameter.
  int reset_sweeps() const { return reset_sweeps_; }

  int num_iterations_to_complete() const override { return min_sweeps_;}
  void set_num_iterations_to_complete(const int sweeps) override;
  int num_iterations() const override { return num_sweeps_; }
  const LnProbability& ln_prob() const override {
    return ln_prob_; }
  void resize(const Histogram& histogram) override;
  std::string write() const override;
  std::string write_per_bin(const int bin) const override;
  std::string write_per_bin_header() const override;
  void set_ln_prob(const LnProbability& ln_prob) override;
  void infrequent_update(const Macrostate& macro) override;
  bool is_equal(const TransitionMatrix& transition_matrix,
    const double tolerance) const;

  /// Return the collection matrix.
  const TripleBandedCollectionMatrix& collection() const {
    return collection_; }

  // HWH hackish interface. See CollectionMatrixSplice::adjust_bounds.
  void set_cm(const int macro, const Bias& bias) override;
  const TripleBandedCollectionMatrix& cm() const override {
    return collection_; }
  const int visits(const int macro) const override { return visits_[macro]; }

  std::shared_ptr<Bias> create(std::istream& istr) const override;
  std::shared_ptr<Bias> create(argtype * args) const override {
    return std::make_shared<TransitionMatrix>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TransitionMatrix(std::istream& istr);
  explicit TransitionMatrix(const Bias& bias);
  virtual ~TransitionMatrix() {}

 private:
  TripleBandedCollectionMatrix collection_;
  LnProbability ln_prob_;
  std::vector<int> visits_;
  int min_visits_ = 0;
  int num_sweeps_ = 0;
  int min_sweeps_ = 0;
  int reset_sweeps_ = -1;
  int average_visits_ = 0;
  int new_sweep_ = 0;
};

inline std::shared_ptr<TransitionMatrix> MakeTransitionMatrix(
    argtype args = argtype()) {
  return std::make_shared<TransitionMatrix>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_TRANSITION_MATRIX_H_
