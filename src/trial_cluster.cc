#include "./trial_cluster.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

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

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

