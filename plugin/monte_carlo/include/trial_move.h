
#ifndef FEASST_MONTE_CARLO_TRIAL_MOVE_H_
#define FEASST_MONTE_CARLO_TRIAL_MOVE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "monte_carlo/include/perturb_move.h"

namespace feasst {

/// Attempt to rigidly move a selection in a Trial in one stage.
class TrialMove : public Trial {
 public:
  TrialMove(
    std::shared_ptr<TrialSelect> select,
    std::shared_ptr<PerturbMove> perturb,
    const argtype& args = argtype()) : Trial(args) {
    add_stage(select, perturb, args);
    set(std::make_shared<TrialComputeMove>());
    class_name_ = "TrialMove";
  }
  explicit TrialMove(std::istream& istr);
  virtual ~TrialMove() {}

 protected:
  void serialize_trial_move_(std::ostream& ostr) const;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_MOVE_H_
