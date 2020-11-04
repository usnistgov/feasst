
#ifndef FEASST_CLUSTER_TRIAL_AVB2_H_
#define FEASST_CLUSTER_TRIAL_AVB2_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"
#include "system/include/neighbor_criteria.h"

namespace feasst {

/// Attempt AVB2 in->out or out->in
std::shared_ptr<Trial> MakeTrialAVB2Half(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    /**
      In addition to the usual Trial arguments:
      - out_to_in: if true, use out->in. Otherwise, in->out.
     */
    const argtype &args = argtype());

/// Attempt AVB2 in->out and out->in with equal probability.
std::shared_ptr<TrialFactory> MakeTrialAVB2(
  std::shared_ptr<NeighborCriteria> neighbor_criteria,
  const argtype &args = argtype());

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_AVB2_H_
