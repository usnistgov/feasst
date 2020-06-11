#include "utils/include/serialize.h"
#include "cluster/include/trial_add_avb.h"
#include "cluster/include/trial_remove_avb.h"
#include "cluster/include/trial_transfer_avb.h"

namespace feasst {

TrialTransferAVB::TrialTransferAVB(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args) : TrialFactory() {
  Arguments args_(args);
  args_.dont_check();
  set_weight(args_.key("weight").dflt("1.").dble());
  add(MakeTrialAddAVB(neighbor_criteria, args));
  add(MakeTrialRemoveAVB(neighbor_criteria, args));
}

}  // namespace feasst
