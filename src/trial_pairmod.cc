#include "./trial_pairmod.h"

namespace feasst {

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

}  // namespace feasst

