
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
std::shared_ptr<Trial> MakeTrialAVB4(const argtype &args = argtype());

// Process AVB4 args, which can also be used in TrialGrow
void gen_avb4_args_(const argtype& args, argtype * args_sel, argtype * args_mv);

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_AVB4_H_
