#include "utils/include/serialize.h"
#include "ewald/include/trial_add_multiple.h"
#include "ewald/include/trial_remove_multiple.h"
#include "ewald/include/trial_transfer_multiple.h"

namespace feasst {

TrialTransferMultiple::TrialTransferMultiple(const argtype& args)
  : TrialFactory(args) {
  add(MakeTrialAddMultiple(args));
  add(MakeTrialRemoveMultiple(args));
}

}  // namespace feasst
