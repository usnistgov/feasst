
#ifndef FEASST_CONFINEMENT_TRIAL_ANYWHERE__H_
#define FEASST_CONFINEMENT_TRIAL_ANYWHERE__H_

#include <memory>
#include "monte_carlo/include/trial.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Attempt to rigidly move anywhere in the box with any orientation.
 */
std::shared_ptr<Trial> MakeTrialAnywhere(argtype args = argtype());

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_TRIAL_ANYWHERE__H_
