
#ifndef FEASST_CLUSTER_TRIAL_TRANSFER_AVB_H_
#define FEASST_CLUSTER_TRIAL_TRANSFER_AVB_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/// Attempt to add a particle with AVB as described in ComputeAddAVB.
std::shared_ptr<Trial> MakeTrialAddAVB(const argtype &args = argtype());

/// Attempt to remove a particle with AVB as described in ComputeRemoveAVB.
std::shared_ptr<Trial> MakeTrialRemoveAVB(const argtype &args = argtype());

/// Attempt TrialAddAVB or TrialRemoveAVB with equal probability.
std::shared_ptr<TrialFactory> MakeTrialTransferAVB(
  const argtype &args = argtype());

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_TRANSFER_AVB_H_
