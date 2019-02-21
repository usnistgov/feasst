
#ifndef FEASST_CORE_TRIAL_H_
#define FEASST_CORE_TRIAL_H_

#include "core/include/perturb.h"
#include "core/include/system.h"
#include "core/include/criteria.h"

namespace feasst {

class Trial {
 public:
  Trial() {
    reset_stats();
    set_weight();
    set_tunable_acceptance();
    set_tunable_percent_change();
  }

  virtual void attempt(Criteria* criteria, System * system) = 0;

  void set_weight(const double weight = 1.) { weight_ = weight; }
  double weight() const { return weight_; }

  void increment_num_attempts() { ++num_attempts_; }
  int64_t num_attempts() const { return num_attempts_; }
  void increment_num_success() { ++num_success_; }
  int64_t num_success() const { return num_success_; }

  /// Set the group index for selection of particle to translate.
  /// By default, a group index of 0 draws from the entire configuration.
  void set_group_index(const int index = 0) { group_index_ = index; }

  /// Return the group index.
  int group_index() const { return group_index_; }

  /// Record a successful attempt.
  void record_success() {
    increment_num_success();
  }

  /// Return the ration of the number of successful attempts and total attempts.
  double acceptance() const {
    return static_cast<double>(num_success_)/static_cast<double>(num_attempts_);
  }

  /// Reset trial statistics.
  virtual void reset_stats() {
    num_attempts_ = 0;
    num_success_ = 0;
  }

  /** @name Tunable
    Tunable parameters to obtain desired acceptance. */
  //@{

  /// Return true if the Trial is tunable.
  bool is_tunable() const { return is_tunable_; }

  /// Set the value of the parameter which is tuned (also enables tunning).
  void set_tunable_param(const double param) {
    is_tunable_ = true;
    tunable_param_ = param; }
  void set_tunable_param_bounded(const double tunable_param) {
    bool is_in_bounds = true;
    if (tunable_param < tunable_param_min_ ||
        tunable_param > tunable_param_max_) {
      is_in_bounds = false;
    }

    if (is_tunable_ && is_in_bounds) {
      set_tunable_param(tunable_param);
    }
  }

  /// Return the above parameter.
  double tunable_param() const { return tunable_param_; }

  /// Set the maximum value of the tunable parameter.
  void set_tunable_param_max(const double max) { tunable_param_max_ = max; }

  /// Return the above parameter.
  double tunable_param_max() const { return tunable_param_max_; }

  /// Set the minimum value of the tunable parameter.
  void set_tunable_param_min(const double min) { tunable_param_min_ = min; }

  /// Return the above parameter.
  double tunable_param_min() const { return tunable_param_min_; }

  /// Set the value of the desired Trial acceptance.
  void set_tunable_acceptance(const double acceptance = 0.25) {
    tunable_acceptance_ = acceptance; }

  /// Return the above parameter.
  double tunable_acceptance() const { return tunable_acceptance_; }

  /// Set the percentage change amount for each tunning step.
  /// Note that a positive percentage means that an increase in the tunable
  /// parameter is expected to result in a decrease in the acceptance.
  /// Input a negative percentage for the reverse expectation.
  virtual void set_tunable_percent_change(const double percent = 0.05) {
    ASSERT(std::abs(percent) < 1., "percent as decimal should be less than 1");
    tunable_percent_change_ = percent;
  }

  /// Return the above parameter
  double tunable_percent_change() const { return tunable_percent_change_; }

  /// Tune the parameter and reset the acceptance statistics.
  virtual void tune() {
    if (is_tunable()) {
      double param = tunable_param();
      if (acceptance() < tunable_acceptance()) {
        param *= 1. - tunable_percent_change_;
      } else {
        param *= 1. + tunable_percent_change_;
      }
      set_tunable_param_bounded(param);
      reset_stats();
    }
  }

  //@}

  /// Return the header description for the status of the trial (e.g., acceptance, etc).
  virtual std::string status_header() const {
    std::stringstream ss;
    ss << "acceptance";
    if (is_tunable()) {
      ss << " tune_param";
    }
    return ss.str();
  }

  /// Return the status of the trial (e.g., acceptance, etc).
  virtual std::string status() const {
    std::stringstream ss;
    ss << acceptance();
    if (is_tunable()) {
      ss << " " << tunable_param();
    }
    return ss.str();
  }

  /// Initialize some variables before each attempt.
  void before_attempt(Criteria* criteria, System * system, Perturb * perturb) {
    perturb->before_attempt();
    criteria->before_attempt(system);
    increment_num_attempts();
  }

  /// Decide to accept or reject the trial, and implement this decision.
  void accept_or_reject(const AcceptanceCriteria& accept_criteria,
      Perturb * perturb,
      Criteria * criteria) {
    if (criteria->is_accepted(accept_criteria)) {
      DEBUG("accepted");
      record_success();
    } else {
      DEBUG("rejected");
      perturb->revert();
    }
  }

  virtual ~Trial() {}

 private:
  double weight_ = 1.;
  int64_t num_attempts_ = 0, num_success_ = 0;
  int group_index_;

  bool is_tunable_ = false;
  double tunable_param_, tunable_acceptance_, tunable_percent_change_;
  double tunable_param_max_, tunable_param_min_;
};

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_H_
