/**
 * \file
 *
 * \brief acceptance criteria for monte carlo trials
 *
 */

#ifndef CRITERIAMETROPOLIS_H_
#define CRITERIAMETROPOLIS_H_

#include "./criteria.h"

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

#endif  // CRITERIAMETROPOLIS_H_

