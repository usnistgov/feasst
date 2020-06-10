
#ifndef FEASST_CLUSTER_TRIAL_AVB2_H_
#define FEASST_CLUSTER_TRIAL_AVB2_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"
#include "system/include/neighbor_criteria.h"

namespace feasst {

/// Attempt AVB2 with both in->out and out->in with equal probability.
class TrialAVB2 : public TrialFactory {
 public:
  explicit TrialAVB2(std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args = argtype());
};

inline std::shared_ptr<TrialAVB2> MakeTrialAVB2(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args = argtype()) {
  return std::make_shared<TrialAVB2>(neighbor_criteria, args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_AVB2_H_
