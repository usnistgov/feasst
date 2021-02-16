
#ifndef FEASST_MONTE_CARLO_TRIAL_REMOVE_MULTIPLE_H_
#define FEASST_MONTE_CARLO_TRIAL_REMOVE_MULTIPLE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {

/**
  Attempt to remove multiple particles.
  Typically requires the use of a reference index.

  args:
  - particle_type[i]: the i-th type of particle to add.
    The "[i]" is to be substituted for an integer 0, 1, 2, ...
 */
std::shared_ptr<Trial> MakeTrialRemoveMultiple(
    argtype args = argtype());

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_REMOVE_MULTIPLE_H_
