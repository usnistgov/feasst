#include "./trial_gca.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

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

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_


