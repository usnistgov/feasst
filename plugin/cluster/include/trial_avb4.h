
#ifndef FEASST_CLUSTER_TRIAL_AVB4_H_
#define FEASST_CLUSTER_TRIAL_AVB4_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {

/**
  Only implemented for single site particles.
  See ComputeAVB4.

  args:
    - neighbor_index: NeighborCriteria index contained in System (default: 0).
 */
std::shared_ptr<Trial> MakeTrialAVB4(argtype args = argtype());

void gen_avb4_args_(argtype * args);

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_AVB4_H_
