#ifndef FEASST_MONTE_CARLO_TRIALS_H_
#define FEASST_MONTE_CARLO_TRIALS_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/// Attempt a rigid translation of a random particle.
std::shared_ptr<Trial> MakeTrialTranslate(argtype args = argtype());

/// Attempt a rigid rotation of a random particle.
std::shared_ptr<Trial> MakeTrialRotate(argtype args = argtype());

/// Attempt to add a particle.
std::shared_ptr<Trial> MakeTrialAdd(argtype args = argtype());

/// Attempt to remove a particle.
std::shared_ptr<Trial> MakeTrialRemove(argtype args = argtype());

/// Attempt TrialAdd or TrialRemove with equal probability.
std::shared_ptr<TrialFactory> MakeTrialTransfer(
  argtype args = argtype());

/// Attempt to change the volume.
std::shared_ptr<Trial> MakeTrialVolume(argtype args = argtype());

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIALS_H_
