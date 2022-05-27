#include "utils/include/serialize.h"
#include "charge/include/trial_add_multiple.h"
#include "charge/include/trial_remove_multiple.h"
#include "charge/include/trial_transfer_multiple.h"

namespace feasst {

class MapTrialTransferMultiple {
 public:
  MapTrialTransferMultiple() {
    auto obj = MakeTrialTransferMultiple();
    obj->deserialize_map()["TrialTransferMultiple"] = obj;
  }
};

static MapTrialTransferMultiple mapper_ = MapTrialTransferMultiple();

TrialTransferMultiple::TrialTransferMultiple(argtype * args) : TrialFactoryNamed() {
  class_name_ = "TrialTransferMultiple";
  argtype orig_args = *args;
  str("shift", args, ""); // remove shift
  auto trial_add = std::make_shared<TrialAddMultiple>(args);
  trial_add->set_weight(trial_add->weight()/2.);
  add(trial_add);
  auto trial_remove = MakeTrialRemoveMultiple(orig_args);
  trial_remove->set_weight(trial_remove->weight()/2.);
  add(trial_remove);
}
TrialTransferMultiple::TrialTransferMultiple(argtype args) : TrialTransferMultiple(&args) {
  FEASST_CHECK_ALL_USED(args);
}

}  // namespace feasst
