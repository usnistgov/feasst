
#ifndef FEASST_CONFINEMENT_TRIAL_ANYWHERE__H_
#define FEASST_CONFINEMENT_TRIAL_ANYWHERE__H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_move.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_anywhere.h"

namespace feasst {

/**
  Attempt to rigidly move anywhere in the box with any orientation.
 */
inline std::shared_ptr<Trial> MakeTrialAnywhere(
    argtype args = argtype()) {
  auto trial = MakeTrialMove(std::make_shared<TrialSelectParticle>(&args),
    std::make_shared<PerturbAnywhere>(),
    "TrialAnywhere",
    &args);
  check_all_used(args);
  return trial;
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_TRIAL_ANYWHERE__H_
