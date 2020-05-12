
#ifndef FEASST_STEPPERS_UTILS_H_
#define FEASST_STEPPERS_UTILS_H_

#include "utils/include/arguments.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/tuner.h"

namespace feasst {

class MonteCarlo;

/**
  Add common steppers to Monte Carlo.
  Includes Log, Movie, CheckEnergy and Tuner.

  args:
  - steps_per: default 1e6
  - file_append: append this to the output files.
  - tolerance: energy check tolerance (default: 1e-10)
 */
void add_common_steppers(MonteCarlo * monte_carlo,
  const argtype& args = argtype());

}  // namespace feasst

#endif  // FEASST_STEPPERS_UTILS_H_
