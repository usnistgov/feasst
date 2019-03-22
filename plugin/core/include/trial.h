
#ifndef FEASST_CORE_TRIAL_H_
#define FEASST_CORE_TRIAL_H_

#include <vector>
#include <numeric>
#include "core/include/perturb.h"
#include "core/include/system.h"
#include "core/include/criteria.h"
#include "core/include/arguments.h"

namespace feasst {

class Trial {
 public:
  Trial(
    /**
      weight: proportional to the likelihood the trial will be attempted
     */
    const argtype &args = argtype()) {
    // default
    reset_stats();
    set_weight();
    set_group_index();

    // parse
    args_.init(args);
    if (args_.key("weight").used()) {
      set_weight(args_.dble());
    }
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

  /// Return the ratio of the number of successful attempts and total attempts.
  double acceptance() const {
    return static_cast<double>(num_success_)/static_cast<double>(num_attempts_);
  }

  /// Reset trial statistics.
  virtual void reset_stats() {
    num_attempts_ = 0;
    num_success_ = 0;
  }

  /// Return the header description for the status of the trial (e.g., acceptance, etc).
  virtual std::string status_header() const {
    std::stringstream ss;
    ss << "acceptance";
    return ss.str();
  }

  /// Return the status of the trial (e.g., acceptance, etc).
  virtual std::string status() const {
    std::stringstream ss;
    ss << acceptance();
    return ss.str();
  }

  /// Initialize some variables before each attempt.
  virtual void before_attempt(Criteria* criteria,
      System * system,
      Perturb * perturb,
      AcceptanceCriteria * accept_criteria) {
    perturb->before_attempt();
    criteria->before_attempt(system);
    increment_num_attempts();
    accept_criteria->ln_metropolis_prob = 0.;
    accept_criteria->energy_new = 0.;
    accept_criteria->energy_new_select = 0.;
    accept_criteria->force_rejection = 0;
    accept_criteria->system = NULL;
    accept_criteria->accepted = -1;
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

  virtual void tune() {}

   virtual void precompute(const std::shared_ptr<Criteria> criteria,
     const System& system) {}

  virtual ~Trial() {}

 protected:
  Arguments args_;

 private:
  double weight_ = 1.;
  int64_t num_attempts_ = 0, num_success_ = 0;
  int group_index_;
};

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_H_
