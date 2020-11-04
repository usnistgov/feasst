
#ifndef FEASST_CONFINEMENT_TRIAL_ANYWHERE_NEW_ONLY_H_
#define FEASST_CONFINEMENT_TRIAL_ANYWHERE_NEW_ONLY_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_move.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_anywhere.h"
#include "mayer/include/trial.h"

namespace feasst {

/**
  Attempt to rigidly move anywhere in the box with any orientation.
  Do not compute the energy of the old configuration.
 */
inline std::shared_ptr<Trial> MakeTrialAnywhereNewOnly(
    const argtype &args = argtype()) {
  auto trial = MakeTrialMove(std::make_shared<TrialSelectParticle>(args),
    std::make_shared<PerturbAnywhere>(),
    "TrialAnywhereNewOnly",
    args);
  trial->set_new_only(true);
  trial->set(std::make_shared<TrialComputeMoveNewOnly>());
  return trial;
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_TRIAL_ANYWHERE_NEW_ONLY_H_
