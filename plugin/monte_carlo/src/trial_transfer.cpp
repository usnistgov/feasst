#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_remove.h"
#include "monte_carlo/include/trial_transfer.h"

namespace feasst {

TrialTransfer::TrialTransfer(const argtype& args) : TrialFactory() {
  Arguments args_(args);
  args_.dont_check();
  set_weight(args_.key("weight").dflt("1.").dble());
  add(MakeTrialAdd(args));
  add(MakeTrialRemove(args));
}

}  // namespace feasst
