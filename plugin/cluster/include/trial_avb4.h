
#ifndef FEASST_CLUSTER_TRIAL_AVB4_H_
#define FEASST_CLUSTER_TRIAL_AVB4_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "system/include/neighbor_criteria.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialAVB4(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args = argtype());

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_AVB4_H_
