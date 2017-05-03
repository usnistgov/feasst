/**
 * \file
 *
 * \brief trial moves for Monte Carlo
 */

#include "./criteria_metropolis.h"
#include "./functions.h"

namespace feasst {

/**
 * Constructor for pair class requires the following
 */
CriteriaMetropolis::CriteriaMetropolis(
  const double beta,   //!< inverse temperature
  const double activ)  //!< activity
    : Criteria(beta, activ) {
  defaultConstruction();
}
/**
 * construct from file
 */
CriteriaMetropolis::CriteriaMetropolis(const char* fileName)
  : Criteria(fileName) {
  defaultConstruction();
}

/**
 * defaults in constructor
 */
void CriteriaMetropolis::defaultConstruction() {
  className_.assign("CriteriaMetropolis");
  verbose_ = 0;
}

CriteriaMetropolis::~CriteriaMetropolis() {
}

/**
 * clone design pattern
 */
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

}  // namespace feasst


