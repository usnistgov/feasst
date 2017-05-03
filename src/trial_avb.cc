#include "./trial_avb.h"

namespace feasst {

TrialAVB::TrialAVB(
  const double pBias,   //!< bias probability
  const double rAbove,   //!< upper limit of bond
  const double rBelow,   //!< lower limit of bond
  const int avbType)  //!< type of avb move
  : Trial(),
    pBias_(pBias),
    avbType_(avbType) {
  initAVB(rAbove, rBelow);
}
TrialAVB::TrialAVB(
  Space *space,
  Pair *pair,
  Criteria *criteria,
  const double pBias,   //!< bias probability
  const double rAbove,   //!< upper limit of bond
  const double rBelow,   //!< lower limit of bond
  const int avbType)   //!< type of avb move
  : Trial(space, pair, criteria),
    pBias_(pBias),
    avbType_(avbType) {
  initAVB(rAbove, rBelow);
}
TrialAVB::TrialAVB(const char* fileName,
                   Space *space,
                   Pair *pair,
                   Criteria *criteria)
  : Trial(space, pair, criteria, fileName) {
  initAVB(rAbove_, rBelow_);
}

}  // namespace feasst

