#ifndef FEASST_MONTE_CARLO_TRIALS_H_
#define FEASST_MONTE_CARLO_TRIALS_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/// Attempt a rigid translation of a random particle.
std::shared_ptr<Trial> MakeTrialTranslate(const argtype &args = argtype());

/// Attempt a rigid rotation of a random particle.
std::shared_ptr<Trial> MakeTrialRotate(const argtype &args = argtype());

/// Attempt to add a particle.
std::shared_ptr<Trial> MakeTrialAdd(const argtype &args = argtype());

/// Attempt to remove a particle.
std::shared_ptr<Trial> MakeTrialRemove(const argtype &args = argtype());

/// Attempt TrialAdd or TrialRemove with equal probability.
std::shared_ptr<TrialFactory> MakeTrialTransfer(
  const argtype &args = argtype());

/// Attempt to change the volume.
std::shared_ptr<Trial> MakeTrialVolume(const argtype &args = argtype());

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIALS_H_
