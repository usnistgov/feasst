
#ifndef FEASST_MONTE_CARLO_UTILS_H_
#define FEASST_MONTE_CARLO_UTILS_H_

#include "utils/include/arguments.h"
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

/// Initialize an add and remove trial simultaneously with the same arguments.
void add_trial_transfer(MonteCarlo * mc, const argtype& args = argtype());

/**
  Prepare a MonteCarlo for a typical LJ simulation.

  args:
  - cubic_box_length: default 8
  - particle: default forcefield/data.lj
  - translate: add a translation trial for all particles (default: true)
  - steps_per: steps per common stepper: (default: 1e4)
  - lrc: use long range corrections (default: true)
 */
void lennard_jones(MonteCarlo * monte_carlo, const argtype& args = argtype());

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_UTILS_H_
