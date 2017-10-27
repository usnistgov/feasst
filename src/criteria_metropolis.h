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

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Metropolis Monte Carlo acceptance criteria.
 */
class CriteriaMetropolis : public Criteria {
 public:
  /// Constructor
  CriteriaMetropolis(const double beta, const double activ);

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

/// Factor method
shared_ptr<CriteriaMetropolis> makeCriteriaMetropolis(const double beta,
                                                      const double activ);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // CRITERIA_METROPOLIS_H_

