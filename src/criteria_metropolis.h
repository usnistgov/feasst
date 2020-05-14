/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef CRITERIA_METROPOLIS_H_
#define CRITERIA_METROPOLIS_H_

#include "./criteria.h"

namespace feasst {

/**
 * Metropolis Monte Carlo acceptance criteria.
 */
class CriteriaMetropolis : public Criteria {
 public:
  // Constructor
  explicit CriteriaMetropolis(const double beta,
    const argtype &args = argtype());

  // Constructor
  // HWH: Depreciate in favor of above
  explicit CriteriaMetropolis(const double beta, const double activ);

  /// Construct by checkpoint file.
  explicit CriteriaMetropolis(const char* fileName);

  ~CriteriaMetropolis() {}
  CriteriaMetropolis* clone() const;
  shared_ptr<CriteriaMetropolis> cloneShrPtr() const;

  /// Return whether to accept (1) or reject (0).
  int accept(const double lnpMet, const double peNew, const char* moveType,
    const int reject) {
    if ( (reject != 1) && (uniformRanNum() < exp(lnpMet)) ) return 1;
    return 0;
    (void) peNew; (void) moveType;  // avoid warning for unused parameters
  }

 protected:
  /// defaults in constructor
  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Criteria> cloneImpl_() const;
};

/// Factory method
shared_ptr<CriteriaMetropolis> makeCriteriaMetropolis(const double beta,
                                                      const double activ);

/// Factory method
shared_ptr<CriteriaMetropolis> makeCriteriaMetropolis(
  const argtype &args = argtype());

}  // namespace feasst

#endif  // CRITERIA_METROPOLIS_H_

