
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

  virtual ~Trial() {}

 private:
  double weight_ = 1.;
  int64_t num_attempts_ = 0;
};

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_H_
