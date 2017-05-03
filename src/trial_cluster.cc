#include "./trial_cluster.h"

namespace feasst {

TrialCluster::TrialCluster(
  const char* transType)   //!< type of transformation
  : Trial(),
    transType_(transType) {
}
TrialCluster::TrialCluster(Space *space,
  Pair *pair,
  Criteria *criteria,
  const char* transType)   //!< type of transformation
  : Trial(space, pair, criteria),
    transType_(transType) {
}
TrialCluster::TrialCluster(const char* fileName,
  Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria, fileName) {
}

}  // namespace feasst

