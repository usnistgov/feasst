/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef TRIAL_DELETE_H_
#define TRIAL_DELETE_H_

#include <memory>
#include <string>
#include "./trial.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class Space;
class Pair;
class Criteria;

/**
 * Attempt to delete particle(s).
 */
class TrialDelete : public Trial {
 public:
  /// Constructor
  TrialDelete(Pair *pair, Criteria *criteria,
              const char* molType);

  /// This constructor is not often used, but its purpose is to initialize trial
  /// for interface before using reconstruct to set object pointers.
  explicit TrialDelete(const char* molType);

  /// Constructors with no molType listed will attempt to delete any particle.
  TrialDelete(Pair *pair, Criteria *criteria);

  /// This constructor is not often used, but its purpose is to initialize trial
  /// for interface before using reconstruct to set object pointers.
  TrialDelete();

  /// Write restart file.
  void writeRestart(const char* fileName);

  /// Construct from restart file.
  TrialDelete(const char* fileName, Pair *pair,
              Criteria *criteria);

  ~TrialDelete() {}
  TrialDelete* clone(Pair *pair, Criteria *criteria) const {
    TrialDelete* t = new TrialDelete(*this);
    t->reconstruct(pair, criteria); return t; }
  shared_ptr<TrialDelete> cloneShrPtr
    (Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialDelete, Trial>
    (cloneImpl(pair, criteria))); }

  // experimental attempt to delete two particles
  int delTwo = 0;

 protected:
  string molType_;   //!< type of molecule to delete
  int molid_;        //!< index of molecule type

  /// attempt deletion
  void attempt1_();

  // default constructor
  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl
    (Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialDelete> t = make_shared<TrialDelete>(*this);
    t->reconstruct(pair, criteria); return t;
  }
};

shared_ptr<TrialDelete> makeTrialDelete(Pair *pair,
  Criteria *criteria, const char* molType);

shared_ptr<TrialDelete> makeTrialDelete(const char* molType);

shared_ptr<TrialDelete> makeTrialDelete(Pair *pair,
  Criteria *criteria);

shared_ptr<TrialDelete> makeTrialDelete();

class MC;

/// Add a "TrialDelete" object to the Monte Carlo object, mc
void deleteTrial(MC *mc, const char* moltype);

/// Add a "TrialDelete" object to the Monte Carlo object, mc
void deleteTrial(shared_ptr<MC> mc, const char* moltype);

/// Add a "TrialDelete" object to the Monte Carlo object, mc
void deleteTrial(MC *mc);

/// Add a "TrialDelete" object to the Monte Carlo object, mc
void deleteTrial(shared_ptr<MC> mc);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // TRIAL_DELETE_H_

