#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_remove.h"
#include "monte_carlo/include/trial_transfer.h"

namespace feasst {

class MapTrialTransfer {
 public:
  MapTrialTransfer() {
    auto obj = MakeTrialTransfer();
    obj->deserialize_map()["TrialTransfer"] = obj;
  }
};

static MapTrialTransfer mapper_ = MapTrialTransfer();

TrialTransfer::TrialTransfer(argtype * args) : TrialFactoryNamed() {
  class_name_ = "TrialTransfer";
  argtype orig_args = *args;
  auto trial_add = std::make_shared<TrialAdd>(args);
  trial_add->set_weight(trial_add->weight()/2.);
  add(trial_add);
  auto trial_remove = MakeTrialRemove(orig_args);
  trial_remove->set_weight(trial_remove->weight()/2.);
  add(trial_remove);
}
TrialTransfer::TrialTransfer(argtype args) : TrialTransfer(&args) {
  FEASST_CHECK_ALL_USED(args);
}

}  // namespace feasst
