/**
 * \file
 *
 * \brief
 */

#include "./trial_pressure.h"

/**
 * Constructor
 */
TrialPressure::TrialPressure(
  const char* variable)    //!< type of transformation
  : Trial(),
    variable_(variable) {
}
TrialPressure::TrialPressure(Space *space,
  Pair *pair,
  Criteria *criteria,
  const char* variable)    //!< type of transformation
  : Trial(space, pair, criteria),
    variable_(variable) {
}
TrialPressure::TrialPressure(const char* fileName,
  Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria, fileName) {
}

