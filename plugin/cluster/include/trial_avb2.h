
#ifndef FEASST_CLUSTER_TRIAL_AVB2_H_
#define FEASST_CLUSTER_TRIAL_AVB2_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Attempt AVB2 in->out or out->in
  Only implemented for single-site particles
  TrialAVB2 below is recommended in most cases to ensure detailed-balance is
  satisfied.
  But this is used in special cases like Prefetch when avoiding TrialFactory.

  args:
    - out_to_in: if true, use out->in. Otherwise, in->out.
    - neighbor_index: NeighborCriteria index contained in System (default: 0).
 */
std::shared_ptr<Trial> MakeTrialAVB2Half(const argtype &args = argtype());

/// Attempt AVB2 in->out and out->in with equal probability.
/// Only implemented for single-site particles
/// See ComputeAVB2 for more information.
std::shared_ptr<TrialFactory> MakeTrialAVB2(const argtype &args = argtype());

// Process AVB2 args, which can also be used in TrialGrow
void gen_avb2_args_(const argtype& args, argtype * args_sel, argtype * args_mv,
  argtype * args_comp);

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_AVB2_H_
