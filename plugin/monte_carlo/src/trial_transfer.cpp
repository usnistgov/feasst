#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_remove.h"
#include "monte_carlo/include/trial_transfer.h"

namespace feasst {

FEASST_MAPPER(TrialTransfer,);

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
  feasst_check_all_used(args);
}

}  // namespace feasst
