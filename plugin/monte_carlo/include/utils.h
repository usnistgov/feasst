
#ifndef FEASST_MONTE_CARLO_UTILS_H_
#define FEASST_MONTE_CARLO_UTILS_H_

#include "utils/include/arguments.h"
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

/// Initialize an add and remove trial simultaneously with the same arguments.
void add_trial_transfer(MonteCarlo * mc, const argtype& args = argtype());

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_UTILS_H_
