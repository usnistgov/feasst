#include "utils/include/serialize.h"
#include "monte_carlo/include/perturb_remove.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_add_avb.h"
#include "cluster/include/compute_add_avb.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/compute_remove_avb.h"
#include "cluster/include/trial_transfer_avb.h"
#include "cluster/include/trial_add_avb.h"
#include "cluster/include/trial_remove_avb.h"

namespace feasst {

class MapTrialTransferAVB {
 public:
  MapTrialTransferAVB() {
    auto obj = MakeTrialTransferAVB();
    obj->deserialize_map()["TrialTransferAVB"] = obj;
  }
};

static MapTrialTransferAVB mapper_trial_transfer_avb__ = MapTrialTransferAVB();

TrialTransferAVB::TrialTransferAVB(argtype * args) : TrialFactoryNamed() {
  class_name_ = "TrialTransferAVB";
  argtype orig_args = *args;
  auto trial_add = MakeTrialAddAVB(orig_args);
  trial_add->set_weight(trial_add->weight()/2.);
  add(trial_add);
  auto trial_remove = std::make_shared<TrialRemoveAVB>(args);
  trial_remove->set_weight(trial_remove->weight()/2.);
  add(trial_remove);
}
TrialTransferAVB::TrialTransferAVB(argtype args) : TrialTransferAVB(&args) {
  FEASST_CHECK_ALL_USED(args);
}

}  // namespace feasst
