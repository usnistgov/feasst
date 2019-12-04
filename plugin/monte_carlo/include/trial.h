
#ifndef FEASST_MONTE_CARLO_TRIAL_H_
#define FEASST_MONTE_CARLO_TRIAL_H_

#include <vector>
#include <string>
#include <memory>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"
#include "monte_carlo/include/perturb.h"

namespace feasst {

/**
  A trial contains a number of TrialStages.
  The Acceptance is computed as the stages are enacted, and then sent to
  Criteria to decide if the trial is accepted or rejected.
 */
class Trial {
 public:
  /**
    args:
    - weight: unnormalized relative probability of selection of this trial
      with respect to all trials (default: 1).
   */
  explicit Trial(const argtype& args = argtype());

  /// Return the unnormalized relative probability of selection of this trial
  /// with respect to all trials.
  double weight() const { return weight_; }

  /// Add a stage which includes selection and perturbation with arguments.
  void add_stage(
    std::shared_ptr<TrialSelect> select,
    std::shared_ptr<Perturb> perturb,
    const argtype& args = argtype());

  /// Set a stage, as can be done just before each attempt.
  void set(const int index, std::shared_ptr<TrialStage> stage);

  /// Number of stages.
  int num_stages() const { return static_cast<int>(stages_.size()); }

  /// Return a stage.
  const TrialStage * stage(const int index) const {
    return stages_[index].get(); }

  // HWH depreciate
  const std::vector<std::shared_ptr<TrialStage> > stages() const {
    return stages_; }

  /// Number of successful attempts.
  int64_t num_success() const { return num_success_; }

  /// Number of attempts.
  int64_t num_attempts() const { return num_attempts_; }

  /// Increment the number of attempts for acceptance.
  void increment_num_attempts() { ++num_attempts_; }

  /// Return the ratio of the number of successful attempts and total attempts.
  double acceptance() const {
    return static_cast<double>(num_success_)/static_cast<double>(num_attempts_);
  }

  /// Reset trial statistics.
  virtual void reset_stats();

  /// Return the header description for the status of the trial
  /// (e.g., acceptance, etc).
  virtual std::string status_header() const;

  /// Return the status of the trial (e.g., acceptance, etc).
  virtual std::string status() const;

  /// Tune parameters.
  virtual void tune();

  /// Set a Mayer-sampling trial.
  void set_mayer(const bool enabled = false) {
    for (auto stage : stages_) stage->set_mayer(enabled); }

  /// Precompute quantities before simulation for optimization.
  virtual void precompute(Criteria * criteria, System * system);

  /// Set the computation of the trial and acceptance.
  void set(std::shared_ptr<TrialCompute> compute) { compute_ = compute; }

  virtual void before_select(Acceptance * acceptance, Criteria * criteria) {}

  /// Revert all stages in reverse order.
  void revert(System * system);

  // Revert changes to system as if trial was never attempted.
  void revert(const int index, const bool accepted, System * system);

  // Finalize all stages in reverse order.
  void finalize(System * system);

  // Require manual finalization of trials (e.g., Pipeline).
  void set_finalize_delayed(const bool delayed = false) {
    is_finalize_delayed_ = delayed; }

  /// Attempt a trial. Return true if accepted.
  virtual bool attempt(Criteria * criteria, System * system, Random * random);

  /* Checks and hacky additions */

  // Return Acceptance, which is a temporary object.
  const Acceptance& accept() const { return acceptance_; }

  // Check if approximately equal to given trial.
  bool is_equal(const Trial& trial) const;

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Trial> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Trial> >& deserialize_map();
  std::shared_ptr<Trial> deserialize(std::istream& istr);
  virtual ~Trial() {}

 protected:
  std::string class_name_ = "Trial";

  void serialize_trial_(std::ostream& ostr) const;
  Trial(std::istream& istr);

  void add_(std::shared_ptr<TrialStage> stage) {
    stages_.push_back(stage);
    refresh_stages_ptr_();
  }

  TrialStage * get_stage_(const int index) { return stages_[index].get(); }

  void increment_num_success_() { ++num_success_; }
  void decrement_num_success_() { --num_success_; }
  void decrement_num_attempts_() { --num_attempts_; }

 private:
  std::vector<std::shared_ptr<TrialStage> > stages_;
  std::shared_ptr<TrialCompute> compute_;
  double weight_ = 1.;
  int64_t num_attempts_ = 0, num_success_ = 0;
  bool is_finalize_delayed_;

  // temporary or duplicate
  Acceptance acceptance_;
  std::vector<TrialStage*> stages_ptr_;

  void refresh_stages_ptr_();
};

inline std::shared_ptr<Trial> MakeTrial(const argtype& args = argtype()) {
  return std::make_shared<Trial>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_H_
