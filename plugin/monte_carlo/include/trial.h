
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

  // Set the weight.
  void set_weight(const double weight) { weight_ = weight; }

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
  const TrialStage& stage(const int index) const {
    return const_cast<TrialStage&>(*stages_[index]); }

  // HWH depreciate
  const std::vector<std::shared_ptr<TrialStage> > stages() const {
    return stages_; }

  /// Number of successful attempts.
  int64_t num_success() const { return data_.int64_1D()[1]; }

  /// Number of attempts.
  int64_t num_attempts() const { return data_.int64_1D()[0]; }

  /// Number of automatic rejections.
  int64_t num_auto_reject() const { return data_.int64_1D()[2]; }

  /// Increment the number of attempts for acceptance.
  void increment_num_attempts() { *num_attempts_() += 1; }
  void increment_num_auto_reject() { *num_auto_reject_() += 1; }

  /// Return the ratio of the number of successful attempts and total attempts.
  double acceptance() const;

  /// Reset trial statistics.
  virtual void reset_stats();

  /// Return the header description for the status of the trial
  /// (e.g., acceptance, etc).
  virtual std::string status_header() const;

  /// Return the status of the trial (e.g., acceptance, etc).
  virtual std::string status() const;

  /// Tune parameters.
  virtual void tune();

  /// Precompute quantities before simulation for optimization.
  virtual void precompute(Criteria * criteria, System * system);

  /// Set the computation of the trial and acceptance.
  void set(std::shared_ptr<TrialCompute> compute) { compute_ = compute; }

  virtual void before_select(Acceptance * acceptance, Criteria * criteria) {}

  /// Revert all stages in reverse order.
  void revert(System * system);

  // Revert changes to system as if trial was never attempted.
  void revert(const int index, const bool accepted, const bool auto_rejected,
              System * system);

  // Finalize all stages in reverse order.
  void finalize(System * system);

  // Require manual finalization of trials (e.g., Pipeline).
  void set_finalize_delayed(const bool delayed = false) {
    is_finalize_delayed_ = delayed; }

  /// Attempt a trial. Return true if accepted.
  virtual bool attempt(Criteria * criteria, System * system, Random * random);

  /// Return the description, as used in Log
  const std::string& description() const { return description_; }

  /// Set the description.
  void set_description(const std::string& description) {
    description_ = description; }

  /* Checks and hacky additions */

  // Return Acceptance, which is a temporary object.
  const Acceptance& accept() const { return acceptance_; }

  // Check if approximately equal to given trial.
  bool is_equal(const Trial& trial) const;

  // prefetch synchronization
  virtual void synchronize_(const Trial& trial) { data_ = trial.data(); }
  const SynchronizeData& data() const { return data_; }

  // Access to factory of Trial objects.
  virtual const std::vector<std::shared_ptr<Trial> >& trials() const;
  virtual const Trial& trial(const int index) const;

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Trial> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Trial> >& deserialize_map();
  std::shared_ptr<Trial> deserialize(std::istream& istr);
  explicit Trial(std::istream& istr);
  virtual ~Trial() {}

 protected:
  std::string class_name_ = "Trial";
  SynchronizeData data_;
  void serialize_trial_(std::ostream& ostr) const;
  void add_(std::shared_ptr<TrialStage> stage) {
    stages_.push_back(stage);
    refresh_stages_ptr_();
  }
  TrialStage * get_stage_(const int index) { return stages_[index].get(); }
  void increment_num_success_() { *num_success_() += 1; }
  void decrement_num_success_() { *num_success_() -= 1; }
  void decrement_num_attempts_() { *num_attempts_() -= 1; }

 private:
  std::vector<std::shared_ptr<TrialStage> > stages_;
  std::shared_ptr<TrialCompute> compute_;
  double weight_;
  std::string description_ = "Trial";
  int64_t * num_attempts_() { return &((*data_.get_int64_1D())[0]); }
  int64_t * num_success_() { return &((*data_.get_int64_1D())[1]); }
  int64_t * num_auto_reject_() { return &((*data_.get_int64_1D())[2]); }
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
