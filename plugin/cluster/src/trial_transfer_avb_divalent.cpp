#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "cluster/include/trial_transfer_avb_divalent.h"
#include "cluster/include/trial_add_avb_divalent.h"
#include "cluster/include/trial_remove_avb_divalent.h"

namespace feasst {

FEASST_MAPPER(TrialTransferAVBDivalent, argtype({{"particle_type", "0"},
  {"particle_type_a", "1"}, {"particle_type_b", "1"}}));

TrialTransferAVBDivalent::TrialTransferAVBDivalent(argtype * args) : TrialFactoryNamed() {
  class_name_ = "TrialTransferAVBDivalent";
  argtype orig_args = *args;
  auto trial_add = MakeTrialAddAVBDivalent(orig_args);
  trial_add->set_weight(trial_add->weight()/2.);
  add(trial_add);
  auto trial_remove = MakeTrialRemoveAVBDivalent(orig_args);
  trial_remove->set_weight(trial_remove->weight()/2.);
  add(trial_remove);
}
TrialTransferAVBDivalent::TrialTransferAVBDivalent(argtype args) : TrialTransferAVBDivalent(&args) {
  // feasst_check_all_used(args);
}
}  // namespace feasst
