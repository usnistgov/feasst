/**
 * \file
 *
 * \brief
 */

#include "./trial_gca.h"

/**
 * Constructor
 */
TrialGCA::TrialGCA() : Trial() {
}
TrialGCA::TrialGCA(Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria) {
}
TrialGCA::TrialGCA(const char* fileName,
  Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria, fileName) {
  maxMoveParam = fstod("maxMoveParam", fileName);
  targAcceptPer = fstod("targAcceptPer", fileName);
}


