
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
inline std::shared_ptr<Trial> MakeTrialMove(
    std::shared_ptr<TrialSelect> select,
    std::shared_ptr<PerturbMove> perturb,
    const std::string& description,
    const argtype& args = argtype()) {
  auto trial = MakeTrial(args);
  trial->set_description(description);
  trial->add_stage(select, perturb, args);
  trial->set(std::make_shared<TrialComputeMove>());
  return trial;
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_MOVE_H_
