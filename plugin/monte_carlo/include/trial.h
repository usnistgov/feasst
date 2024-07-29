
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
  //@{
  /** @name Arguments
   */

  /**
    args:
    - weight: unnormalized relative probability of selection of this trial
      with respect to all trials (default: 1).
    - weight_per_number_fraction: if > 0, the weight is continuously updated to
      (weight_per_number_fraction * number of TrialSelect::particle_type in
       first TrialStage / total number of particles).
      If <= 0, the given weight above is fixed to that value (default: -1).
      For a multicomponent simulation, weight_per_number_fraction with a separate
      trial for each mobile particle_type ensures equal probability of selecting
      any particles of those types type regardless of the number of particles of
      each type, while still allowing for different tunable parameters for each
      particle_type.
    - number_fraction_exclude_type: if >= 0 (default: -1), exclude this particle
      type from the total number of fractions in weight_per_number_fraction.
      If there is a rigid particle type in the grand canonical ensemble, avoid
      changing the relative weight of translation vs insert/delete and breaking
      detailed balance.
   */
  explicit Trial(argtype args = argtype());
  explicit Trial(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return the unnormalized relative probability of selection of this trial
  /// with respect to all trials.
  double weight() const { return data_.dble_1D()[0]; }

  // Set the weight.
  void set_weight(const double weight) { *get_weight_() = weight; }

  /// Return the weight per number of particles of given type.
  double weight_per_number_fraction() const {
    return weight_per_number_fraction_; }

  // Set the weight per number.
  void set_weight_per_number_fraction(const double wpn) {
    weight_per_number_fraction_ = wpn; }

  /// Return the particle type excluded in number fraction weights.
  int number_fraction_exclude_type() const {
    return number_fraction_exclude_type_; }

  // Set the particle type excluded in number fraction weights.
  void set_number_fraction_exclude_type(const double nfet) {
    number_fraction_exclude_type_ = nfet; }

  /// Add a stage which includes selection and perturbation with arguments.
  void add_stage(std::shared_ptr<TrialSelect> select,
                 std::shared_ptr<Perturb> perturb,
                 argtype * stage_args);

  /// Same as above, but without arguments for stage.
  void add_stage(std::shared_ptr<TrialSelect> select,
                 std::shared_ptr<Perturb> perturb);

  /// Same as above, but by copying an existing stage.
  void add_stage(std::shared_ptr<TrialStage> stage);

  /// Set a stage, as can be done just before each attempt.
  void set(const int index, std::shared_ptr<TrialStage> stage);

  /// Number of stages.
  int num_stages() const { return static_cast<int>(stages_.size()); }

  /// Return a stage.
  const TrialStage& stage(const int index) const {
    return const_cast<TrialStage&>(*stages_[index]); }

  // HWH deprecate
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
  /// Return -1 if no attempts were performed.
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

  /// Return TrialCompute.
  const TrialCompute& compute() const {
    return const_cast<TrialCompute&>(*compute_); }

  virtual void before_select(Acceptance * acceptance, Criteria * criteria) {}

  /// Revert all stages in reverse order.
  void revert(System * system, Criteria * criteria);

  // Revert changes to system as if trial was never attempted.
  void revert(const int index, const bool accepted, const bool auto_rejected,
              System * system, Criteria * criteria);

  // Finalize all stages in reverse order.
  void finalize(System * system, Criteria * criteria);

  // Require manual finalization of trials (e.g., Pipeline).
  void set_finalize_delayed(const bool delayed = false) {
    is_finalize_delayed_ = delayed; }

  /// Attempt a trial. Return true if accepted.
  virtual bool attempt(Criteria * criteria, System * system, Random * random);

  //HWH Depreciate description once all Trials are derived classes again.
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

  TrialStage * get_stage_(const int index) { return stages_[index].get(); }

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Trial> create(std::istream& istr) const;
  virtual std::shared_ptr<Trial> create(argtype * args) const;
  std::map<std::string, std::shared_ptr<Trial> >& deserialize_map();
  std::shared_ptr<Trial> deserialize(std::istream& istr);
  std::shared_ptr<Trial> factory(const std::string name, argtype * args);
  explicit Trial(std::istream& istr);
  virtual ~Trial() {}

  //@}
 protected:
  std::string class_name_ = "Trial";
  SynchronizeData data_;
  void serialize_trial_(std::ostream& ostr) const;
  void increment_num_success_() { *num_success_() += 1; }
  void decrement_num_success_() { *num_success_() -= 1; }
  void decrement_num_attempts_() { *num_attempts_() -= 1; }

 private:
  std::vector<std::shared_ptr<TrialStage> > stages_;
  std::shared_ptr<TrialCompute> compute_;
  //double weight_;
  double * get_weight_() { return &((*data_.get_dble_1D())[0]); }
  double weight_per_number_fraction_;
  int number_fraction_exclude_type_;
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

inline std::shared_ptr<Trial> MakeTrial(argtype args = argtype()) {
  return std::make_shared<Trial>(args); }

inline std::shared_ptr<Trial> MakeTrial(argtype * args) {
  return std::make_shared<Trial>(args); }

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_H_
