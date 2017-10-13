/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

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

shared_ptr<CriteriaMetropolis> makeCriteriaMetropolis(const double beta,
  const double activ) {
  return make_shared<CriteriaMetropolis>(beta, activ);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_


