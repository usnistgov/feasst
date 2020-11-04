#include "utils/include/serialize.h"
#include "ewald/include/trial_add_multiple.h"
#include "ewald/include/trial_remove_multiple.h"
#include "ewald/include/trial_transfer_multiple.h"

namespace feasst {

std::shared_ptr<TrialFactory> MakeTrialTransferMultiple(
    const argtype &args) {
  auto factory = std::make_shared<TrialFactory>(args);
  factory->add(MakeTrialAddMultiple(args));
  factory->add(MakeTrialRemoveMultiple(args));
  return factory;
}

}  // namespace feasst
