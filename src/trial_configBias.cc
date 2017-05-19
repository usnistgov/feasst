#include "./trial_configBias.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

TrialConfigBias::TrialConfigBias(const char* molType)
  : Trial(),
    molType_(molType) {
}
TrialConfigBias::TrialConfigBias(Space *space,
  Pair *pair,
  Criteria *criteria,
  const char* molType)   //!< type of molecule to add
  : Trial(space, pair, criteria),
    molType_(molType) {
  // // this line creates a movie
  // pair_->printxyz("asdf", 1);
}
TrialConfigBias::TrialConfigBias(const char* fileName,
  Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria, fileName) {
}

/**
 * clone design pattern
 */
TrialConfigBias* TrialConfigBias::clone
  (Space* space, Pair *pair, Criteria *criteria) const {
  TrialConfigBias* t = new TrialConfigBias(*this);
  t->reconstruct(space, pair, criteria);
  return t;
}
shared_ptr<TrialConfigBias> TrialConfigBias::cloneShrPtr
  (Space* space, Pair* pair, Criteria* criteria) const {
  return(std::static_pointer_cast<TrialConfigBias, Trial>
    (cloneImpl(space, pair, criteria)));
}
shared_ptr<Trial> TrialConfigBias::cloneImpl
  (Space* space, Pair *pair, Criteria *criteria) const {
  shared_ptr<TrialConfigBias> t = make_shared<TrialConfigBias>(*this);
  t->reconstruct(space, pair, criteria);
  return t;
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

