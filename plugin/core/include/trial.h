
#ifndef FEASST_CORE_TRIAL_H_
#define FEASST_CORE_TRIAL_H_

#include "core/include/perturb.h"
#include "core/include/system.h"
#include "core/include/criteria.h"

namespace feasst {

class Trial {
 public:
  Trial() { reset_stats(); }

  virtual void attempt(Criteria* criteria, System * system) = 0;

  void set_weight(const double weight) { weight_ = weight; }
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

  /// Return the header description for the status of the trial (e.g., acceptance, etc).
  virtual std::string status_header() const {
    return std::string("acceptance");
  }

  /// Return the status of the trial (e.g., acceptance, etc).
  virtual std::string status() const {
    std::stringstream ss;
    ss << static_cast<double>(num_success_)/static_cast<double>(num_attempts_);
    return ss.str();
  }

  /// Reset trial statistics.
  virtual void reset_stats() {
    num_attempts_ = 0;
    num_success_ = 0;
  }

  virtual ~Trial() {}

 private:
  double weight_ = 1.;
  int64_t num_attempts_ = 0;
  int64_t num_success_ = 0;
  int group_index_;
};

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_H_
