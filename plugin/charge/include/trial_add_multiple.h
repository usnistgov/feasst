
#ifndef FEASST_CHARGE_TRIAL_ADD_MULTIPLE_H_
#define FEASST_CHARGE_TRIAL_ADD_MULTIPLE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {


// parse the number of particle types.
std::vector<int> ptypes(argtype * args);

/**
  Attempt to add multiple particles.
  Typically requires the use of a reference index.

  args:
  - particle_type[i]: the i-th type of particle to add.
    The "[i]" is to be substituted for an integer 0, 1, 2, ...
 */
std::shared_ptr<Trial> MakeTrialAddMultiple(argtype args = argtype());

}  // namespace feasst

#endif  // FEASST_CHARGE_TRIAL_ADD_MULTIPLE_H_
