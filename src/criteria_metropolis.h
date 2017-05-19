/**
 * \file
 *
 * \brief acceptance criteria for monte carlo trials
 *
 */

#ifndef CRITERIA_METROPOLIS_H_
#define CRITERIA_METROPOLIS_H_

#include "./criteria.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class CriteriaMetropolis : public Criteria {
 public:
  CriteriaMetropolis(const double beta, const double activ);
  explicit CriteriaMetropolis(const char* fileName);
  ~CriteriaMetropolis();
  CriteriaMetropolis* clone() const;
  shared_ptr<CriteriaMetropolis> cloneShrPtr() const;

  /// defaults in constructor
  void defaultConstruction();

  /// acceptance criteria for trial moves
  int accept(const double lnpMet, const double peNew, const char* moveType,
    const int reject) {
    if ( (reject != 1) && (uniformRanNum() < exp(lnpMet)) ) return 1;
    return 0;
    (void) peNew; (void) moveType;  // unused parameters
  }

 protected:
  // clone design pattern
  virtual shared_ptr<Criteria> cloneImpl_() const;
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // CRITERIA_METROPOLIS_H_

