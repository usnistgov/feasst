
#ifndef FEASST_CORE_TRIAL_H_
#define FEASST_CORE_TRIAL_H_

#include "core/include/perturb.h"
#include "core/include/system.h"
#include "core/include/criteria.h"

namespace feasst {

class Trial {
 public:
  virtual void attempt(Criteria* criteria, System * system) = 0;

  void set_weight(const double weight) { weight_ = weight; }
  double weight() const { return weight_; }

  void increment_num_attempts() { ++num_attempts_; }
  int64_t num_attempts() const { return num_attempts_; }

  /// Set the group index for selection of particle to translate.
  /// By default, a group index of 0 draws from the entire configuration.
  void set_group_index(const int index = 0) { group_index_ = index; }

  /// Return the group index.
  int group_index() const { return group_index_; }

  virtual ~Trial() {}

 private:
  double weight_ = 1.;
  int64_t num_attempts_ = 0;
  int group_index_;
};

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_H_
