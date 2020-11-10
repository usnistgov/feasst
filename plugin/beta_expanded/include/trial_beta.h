#ifndef FEASST_BETA_EXPANDED_TRIAL_BETA_H_
#define FEASST_BETA_EXPANDED_TRIAL_BETA_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {

/// Attempt to change the inverse temperature, \f$\beta\f$ by a fixed amount.
std::shared_ptr<Trial> MakeTrialBeta(const argtype &args = argtype());

}  // namespace feasst

#endif  // FEASST_BETA_EXPANDED_TRIAL_BETA_H_
