
#ifndef FEASST_FLAT_HISTOGRAM_FLAT_HISTOGRAM_H_
#define FEASST_FLAT_HISTOGRAM_FLAT_HISTOGRAM_H_

#include <memory>
#include "utils/include/arguments.h"
#include "math/include/accumulator.h"
#include "monte_carlo/include/criteria.h"
#include "flat_histogram/include/macrostate.h"
#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/bias.h"

namespace feasst {

class Random;

/**
  Flat histogram acceptance criteria uses a bias to improve sampling and attempt
  to recover the free energy of the system as a function of the give macrostate.

  The Macrostate must be defined before the bias.
 */
class FlatHistogram : public Criteria {
 public:
  // HWH remove, considering EnsembleAverage?
  // Only used for reweighting
  FlatHistogram();

  /// Construct with a macrostate and a bias
  FlatHistogram(std::shared_ptr<Macrostate> macrostate,
      std::shared_ptr<Bias> bias);

  /**
    This is a flattened constructor which takes arguments for macrostate, the
    macrostate histogram and the bias, as well as the following two arguments.

    args:
    - macrostate: MacrostateNumParticles, MacrostateEnergy, etc
    - bias: WangLandau, TransitionMatrix, WLTM, etc.
   */
  explicit FlatHistogram(argtype args);
  explicit FlatHistogram(argtype * args);

  /// Same as above, but with an added Constraint.
  FlatHistogram(std::shared_ptr<Macrostate> macrostate,
      std::shared_ptr<Bias> bias,
      std::shared_ptr<Constraint> constraint);

  /// Return the macrostate.
  const Macrostate& macrostate() const {
    return const_cast<Macrostate&>(*macrostate_); }

  /// Return the bias.
  const Bias& bias() const { return const_cast<Bias&>(*bias_); }
  void set_bias(std::shared_ptr<Bias> bias) { bias_ = bias; }

  int num_iterations_to_complete() const override {
    return bias_->num_iterations_to_complete(); }
  void set_num_iterations_to_complete(const int num) override {
    bias_->set_num_iterations_to_complete(num); }
  int num_iterations() const override { return bias_->num_iterations(); }
  bool is_complete() const override { return bias_->is_complete(); }

  void before_attempt(const System& system) override;

  bool is_accepted(const Acceptance& acceptance,
    const System& system,
    Random * random) override;

  std::string write() const override;
  int phase() const override { return bias_->phase(); }
  void increment_phase() override { bias_->increment_phase(); }

  /// Return the state. Return -1 if state is not determined.
  int state() const override { return macrostate_current_; }
  int num_states() const override { return macrostate_->histogram().size(); }
  int state_old() const override { return macrostate_old_; }
  int state_new() const override { return macrostate_new_; }

  /// Set the macrostate probability distribution.
  void set_ln_prob(const LnProbability& ln_prob) {
    bias_->set_ln_prob(ln_prob); }

  /// Return the macrostate probability distribution.
  const LnProbability& ln_prob() const { return bias_->ln_prob(); }

  // HWH hackish implementation for prefetch
  // Revert changes from previous trial.
  // HWH rename: delete
  void revert_(const bool accepted, const double ln_prob) override;
  void imitate_trial_rejection_(const double ln_prob,
      const int state_old,
      const int state_new) override {
    bias_->update(state_old, state_new, ln_prob, false); }

  void update() override { bias_->infrequent_update(); }

  bool is_fh_equal(const FlatHistogram& flat_histogram,
    const double tolerance) const;

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<FlatHistogram>(istr); }
  std::shared_ptr<Criteria> create(argtype * args) const override {
    return std::make_shared<FlatHistogram>(args); }
  void serialize(std::ostream& ostr) const override;
  FlatHistogram(std::istream& istr);
  FlatHistogram(const Criteria& criteria);
  ~FlatHistogram() {}

 private:
  std::shared_ptr<Bias> bias_;
  std::shared_ptr<Macrostate> macrostate_;
  int macrostate_old_ = -1;
  int macrostate_new_ = -1;
  int macrostate_current_ = -1;
  bool is_macrostate_set_ = false;

  // temporary
  Acceptance empty_;

  void init_(std::shared_ptr<Macrostate> macrostate,
    std::shared_ptr<Bias> bias);
};

inline std::shared_ptr<FlatHistogram> MakeFlatHistogram() {
  return std::make_shared<FlatHistogram>();
}

inline std::shared_ptr<FlatHistogram> MakeFlatHistogram(
    std::shared_ptr<Macrostate> macrostate,
    std::shared_ptr<Bias> bias) {
  return std::make_shared<FlatHistogram>(macrostate, bias);
}

inline std::shared_ptr<FlatHistogram> MakeFlatHistogram(argtype args) {
  return std::make_shared<FlatHistogram>(args);
}

inline std::shared_ptr<FlatHistogram> MakeFlatHistogram(
    std::shared_ptr<Macrostate> macrostate,
    std::shared_ptr<Bias> bias,
    std::shared_ptr<Constraint> constraint) {
  return std::make_shared<FlatHistogram>(macrostate, bias, constraint);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_FLAT_HISTOGRAM_H_
