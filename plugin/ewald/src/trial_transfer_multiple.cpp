#include "utils/include/serialize.h"
#include "ewald/include/trial_add_multiple.h"
#include "ewald/include/trial_remove_multiple.h"
#include "ewald/include/trial_transfer_multiple.h"

namespace feasst {

std::shared_ptr<TrialFactory> MakeTrialTransferMultiple(argtype args) {
  argtype rm_args = args, add_args = args;
  str("shift", &add_args, ""); // remove shift
  auto factory = std::make_shared<TrialFactory>(&args);
  factory->add(MakeTrialAddMultiple(add_args));
  factory->add(MakeTrialRemoveMultiple(rm_args));
  return factory;
}

}  // namespace feasst
