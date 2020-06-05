#include "utils/include/utils_io.h"
#include "system/include/utils.h"
#include "monte_carlo/include/utils.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_remove.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"

namespace feasst {

void add_trial_transfer(MonteCarlo * mc, const argtype& args) {
  mc->add(MakeTrialAdd(args));
  mc->add(MakeTrialRemove(args));
}

}  // namespace feasst
