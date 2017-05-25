#include "./criteria_metropolis.h"
#include "./functions.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

CriteriaMetropolis::CriteriaMetropolis(const double beta, const double activ)
  : Criteria(beta, activ) {
  defaultConstruction_();
}

CriteriaMetropolis::CriteriaMetropolis(const char* fileName)
  : Criteria(fileName) {
  defaultConstruction_();
}

void CriteriaMetropolis::defaultConstruction_() {
  className_.assign("CriteriaMetropolis");
  verbose_ = 0;
}

CriteriaMetropolis* CriteriaMetropolis::clone() const {
  CriteriaMetropolis* c = new CriteriaMetropolis(*this);
  c->reconstruct();
  return c;
}

shared_ptr<CriteriaMetropolis> CriteriaMetropolis::cloneShrPtr() const {
  return(std::static_pointer_cast<CriteriaMetropolis, Criteria>(cloneImpl_()));
}

shared_ptr<Criteria> CriteriaMetropolis::cloneImpl_() const {
  shared_ptr<CriteriaMetropolis> c = make_shared<CriteriaMetropolis>(*this);
  c->reconstruct();
  return c;
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_


