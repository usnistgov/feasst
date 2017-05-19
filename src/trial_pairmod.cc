#include "./trial_pairmod.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

TrialPairMod::TrialPairMod() : Trial() {
}
TrialPairMod::TrialPairMod(Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria) {
}
TrialPairMod::TrialPairMod(const char* fileName,
  Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria, fileName) {
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

