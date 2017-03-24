/**
 * \file
 *
 * \brief
 */

#include "./trial_grow.h"

/**
 * Constructor
 */
TrialGrow::TrialGrow(
  const char* molType,   //!< type of molecule to grow
  const int nStages)     //!< number of stages in growth
  : Trial(),
    molType_(molType) {
}
TrialGrow::TrialGrow(Space *space,
  Pair *pair,
  Criteria *criteria,
  const char* molType,   //!< type of molecule to grow
  const int nStages)     //!< number of stages in growth
  : Trial(space, pair, criteria),
    molType_(molType) {
}
TrialGrow::TrialGrow(const char* fileName,
  Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria, fileName) {
}

