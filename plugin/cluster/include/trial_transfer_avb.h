
#ifndef FEASST_CLUSTER_TRIAL_TRANSFER_AVB_H_
#define FEASST_CLUSTER_TRIAL_TRANSFER_AVB_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"
#include "system/include/neighbor_criteria.h"

namespace feasst {

/// Attempt TrialAddAVB or TrialRemoveAVB with equal probability.
class TrialTransferAVB : public TrialFactory {
 public:
  explicit TrialTransferAVB(std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args = argtype());
};

inline std::shared_ptr<TrialTransferAVB> MakeTrialTransferAVB(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args = argtype()) {
  return std::make_shared<TrialTransferAVB>(neighbor_criteria, args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_TRANSFER_AVB_H_
