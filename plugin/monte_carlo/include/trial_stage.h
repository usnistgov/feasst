
#ifndef FEASST_MONTE_CARLO_TRIAL_STAGE_H_
#define FEASST_MONTE_CARLO_TRIAL_STAGE_H_

#include <string>
#include <memory>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/perturb.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/rosenbluth.h"

namespace feasst {

/**
  A stage contains both a selection and perturbation.
  Random realizations of the perturbation are performed in a number of steps.
  A reference potential may be used during these steps, following the dual-cut
  configurational bias approach: http://doi.org/10.1080/002689798167881
 */
class TrialStage {
 public:
  /**
    args:
    - num_steps: number of rosenbluth steps (default: 1).
    - reference_index: index of reference potential.
      Otherwise, if full potential is desired, set to -1 (default: -1).
   */
  explicit TrialStage(const argtype& args = argtype());

  /// Return the index of the reference potential.
  int reference() const { return reference_; }

  /// Return true if the trial is utilizing Mayer sampling.
  bool is_mayer() const { return is_mayer_; }

  /// Set the above.
  void set_mayer(const bool enabled = false) { is_mayer_ = enabled; }

  /// Return the Rosenbluth.
  const Rosenbluth& rosenbluth() const { return rosenbluth_; }

  /// Set the selection.
  void set(std::shared_ptr<TrialSelect> select) { select_ = select; }

  /// Return the above.
  const TrialSelect * trial_select() const { return select_.get(); }

  /// Initialization before any stage attempt.
  void precompute(System * system);

  /// Initializations before each stage attempt.
  void before_select();

  /// Perform the selection and update the acceptance.
  /// Return false if selection fails.
  /// Otherwise, set sites involved in stage as unphysical and return true.
  bool select(System * system, Acceptance * acceptance, Random * random);

  /// Set the perturbation.
  void set(std::shared_ptr<Perturb> perturb) { perturb_ = perturb; }

  /// Return the above.
  const Perturb * perturb() const { return perturb_.get(); }

  /// Set mobile selection physical.
  void set_mobile_physical(const bool physical, System * system);

  /// Attempt all steps in a stage.
  /// Consider reference potentials and compute Rosenbluth factors.
  /// Set sites involves in stage as physical.
  void attempt(
    System * system,
    Acceptance * acceptance,
    Criteria * criteria,
    /// Set to 1 for "old" system and "0" for new.
    const int old,
    Random * random);

  /// Call between multiple attempts (e.g., old vs new)
  /// Set mobile sites unphysical.
  void mid_stage(System * system);

  /// Return true if constraints are satisfied.
  bool are_constraints_satisfied(const System& system) const;

  /// Revert the attempt.
  void revert(System * system) { perturb_->revert(system); }

  /// Finalize the attempt.
  void finalize(System * system) { perturb_->finalize(system); }

  /// Tune parameters.
  void tune(const double acceptance) { perturb_->tune(acceptance); }

  /// Print status header.
  std::string status_header() const { return perturb_->status_header(); }

  /// Print status.
  std::string status() const { return perturb_->status(); }

  // HWH avoid using this
  TrialSelect * get_trial_select() { return select_.get(); }

  /// Serialize.
  void serialize(std::ostream& ostr) const;

  /// Deserialize.
  explicit TrialStage(std::istream& istr);

  ~TrialStage() {}

 private:
  int reference_ = -1;
  std::shared_ptr<Perturb> perturb_;
  std::shared_ptr<TrialSelect> select_;
  Rosenbluth rosenbluth_;
  bool is_mayer_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_STAGE_H_
