
#ifndef FEASST_FLAT_HISTOGRAM_FLAT_HISTOGRAM_H_
#define FEASST_FLAT_HISTOGRAM_FLAT_HISTOGRAM_H_

#include <memory>
#include "math/include/accumulator.h"
#include "monte_carlo/include/criteria.h"
#include "flat_histogram/include/macrostate.h"
#include "flat_histogram/include/bias.h"

namespace feasst {

/**
  The macrostate must be defined before the bias.
  Use MacrostateAccumulator to compute custom per-macrostate quantities.
 */
class FlatHistogram : public Criteria {
 public:
  FlatHistogram(const argtype &args = argtype()) : Criteria(args) {
    class_name_ = "FlatHistogram";
  }

  FlatHistogram(std::shared_ptr<Macrostate> macrostate,
      std::shared_ptr<Bias> bias,
      const argtype &args = argtype()) : FlatHistogram(args) {
//    class_name_ = "FlatHistogram";
    set(macrostate);
    set(bias);
  }

  /// Same as above, but with an added Constraint.
  FlatHistogram(std::shared_ptr<Macrostate> macrostate,
      std::shared_ptr<Bias> bias,
      std::shared_ptr<Constraint> constraint,
      const argtype &args = argtype()) : FlatHistogram(macrostate, bias, args) {
    add(constraint);
  }

  void before_attempt(const System* system) override;

  bool is_accepted(const Acceptance& acceptance,
    const System * system,
    const double uniform_random) override;

  std::string write() const override;
  bool is_complete() const override { return bias_->is_complete(); }

  /// Set the macrostate which is subject to the bias.
  void set(const std::shared_ptr<Macrostate> macrostate) {
    macrostate_ = macrostate;
    is_macrostate_set_ = true;
  }

  /// Return the macrostate.
  const Macrostate * macrostate() const { return macrostate_.get(); }

  /// Set the bias for the flat histogram method.
  void set(const std::shared_ptr<Bias> bias);

  /// Return the bias.
  const Bias * bias() const { return bias_.get(); }

  /// Return the state. Return -1 if state is not determined.
  int state() const override { return macrostate_current_; }
  int num_states() const override { return macrostate_->histogram().size(); }
  int state_old() const override { return macrostate_old_; }
  int state_new() const override { return macrostate_new_; }

  /// Revert changes from previous trial.
  // HWH rename: delete
  void revert(const bool accepted, const double ln_prob) override;

  // HWH consider moving the below functions to a new class or util

  /// Return the pressure.
  /// Return a reweighted macrostate distribution.
  LnProbability reweight(
      /**
        Change in conjugate variable. For example, if reweighting in number
        of particles using the grand canonical ensemble, the change in the
        thermodynamic conjugate is \f$\Delta(\beta\mu\)f$.
        If reweighting in potential energy in the microcanonical, the change
        in the thermodynamic conjugate is \f$\Delta(-\beta)\f$.
       */
      const double delta_conjugate);

  /// Set the macrostate probability distribution.
  void set_ln_prob(const LnProbability& ln_prob) {
    bias_->set_ln_prob(ln_prob); }

  // HWH add reference to Shen/Errington
  double pressure(const double volume,
      /// Select phase by order of macrostate.
      /// Assumes default method of dividing phase boundary.
      const int phase = 0) const;

  /// Return the ensemble averaged property from a list of properties averaged
  /// at each macrostate.
  double average(const LnProbability& ln_prob,
                 // HWH make macrostate_averages work with bin averages.
                 const std::vector<double>& macrostate_averages,
                 /// Select phase by order of macrostate.
                 /// Assumes default method of dividing phase boundary.
                 const int phase = 0) const;

  double average(const std::vector<double>& macrostate_averages,
                 /// Select phase by order of macrostate.
                 /// Assumes default method of dividing phase boundary.
                 const int phase = 0) const {
  return average(ln_prob_(), macrostate_averages, phase); }

  /// Return the average macrostate
  double average_macrostate(const LnProbability& ln_prob,
      /// Select phase by order of macrostate.
      /// Assumes default method of dividing phase boundary.
      const int phase = 0) const;

  double average_macrostate(
      /// Select phase by order of macrostate.
      /// Assumes default method of dividing phase boundary.
      const int phase = 0) const {
    return average_macrostate(ln_prob_(), phase);
  }

  void imitate_trial_rejection(const double ln_prob,
      const int state_old,
      const int state_new) override {
    bias_->update(state_old, state_new, ln_prob, false);
  }

  void update() override { bias_->infrequent_update(); }

  bool is_fh_equal(const FlatHistogram* flat_histogram,
    const double tolerance) const;

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<FlatHistogram>(istr); }
  void serialize(std::ostream& ostr) const override;
  FlatHistogram(std::istream& istr);
  ~FlatHistogram() {}

 private:
  std::shared_ptr<Bias> bias_;
  std::shared_ptr<Macrostate> macrostate_;
  int macrostate_old_ = -1;
  int macrostate_new_ = -1;
  int macrostate_current_ = -1;
  bool is_macrostate_set_ = false;

  const LnProbability& ln_prob_() const { return bias_->ln_prob(); }

  /// Determine min and max indices for a given phase
  /// Return -1 if no phase boundary.
  void phase_boundary_(const LnProbability& ln_prob,
      const int phase, int * min, int * max) const;

  /// Same as above but with the ln_prob_ contained in this class.
  void phase_boundary_(const int phase, int * min, int * max) const {
    phase_boundary_(ln_prob_(), phase, min, max);
  }

  // temporary
  Acceptance empty_;
};

inline std::shared_ptr<FlatHistogram> MakeFlatHistogram(
    const argtype &args = argtype()) {
  return std::make_shared<FlatHistogram>(args);
}

inline std::shared_ptr<FlatHistogram> MakeFlatHistogram(
    std::shared_ptr<Macrostate> macrostate,
    std::shared_ptr<Bias> bias,
    const argtype &args = argtype()) {
  return std::make_shared<FlatHistogram>(macrostate, bias, args);
}

inline std::shared_ptr<FlatHistogram> MakeFlatHistogram(
    std::shared_ptr<Macrostate> macrostate,
    std::shared_ptr<Bias> bias,
    std::shared_ptr<Constraint> constraint,
    const argtype &args = argtype()) {
  return std::make_shared<FlatHistogram>(macrostate, bias, constraint, args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_FLAT_HISTOGRAM_H_
