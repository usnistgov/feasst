
#ifndef FEASST_MONTE_CARLO_TRIAL_TRANSFER_MULTIPLE_H_
#define FEASST_MONTE_CARLO_TRIAL_TRANSFER_MULTIPLE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/// Attempt TrialAddMultiple or TrialRemoveMultiple and split the trial weights
/// equally.
std::shared_ptr<TrialFactory> MakeTrialTransferMultiple(
    const argtype &args = argtype());

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_TRANSFER_MULTIPLE_H_
