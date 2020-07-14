#include "utils/include/serialize.h"
#include "cluster/include/trial_avb2_half.h"
#include "cluster/include/trial_avb2.h"

namespace feasst {

TrialAVB2::TrialAVB2(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args) : TrialFactory(args) {
  argtype out2in_args(args);
  out2in_args.insert({"out_to_in", "true"});
  add(MakeTrialAVB2Half(neighbor_criteria, out2in_args));

  argtype in2out_args(args);
  in2out_args.insert({"out_to_in", "false"});
  add(MakeTrialAVB2Half(neighbor_criteria, in2out_args));
}

}  // namespace feasst
