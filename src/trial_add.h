/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef TRIAL_ADD_H_
#define TRIAL_ADD_H_

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
 * Attempt to add particle(s).
 */
class TrialAdd : public Trial {
 public:
  /// Construction
  /// @param molType description of molecule to add which matches addMol.
  TrialAdd(Pair *pair, Criteria *criteria, const char* molType);

  /// This constructor is not often used, but its purpose is to initialize trial
  /// for interface before using reconstruct to set object pointers.
  explicit TrialAdd(const char* molType);

  /// Write restart file.
  void writeRestart(const char* fileName);

  /// Construct from restart file.
  TrialAdd(const char* fileName, Pair *pair, Criteria *criteria);

  // Overloaded from base class for status of specific trials.
  string printStat(const bool header = false);

  /// read only access to protected variables
  string molType() const { return molType_; }

  ~TrialAdd() {}
  TrialAdd* clone(Pair *pair, Criteria *criteria) const;
  shared_ptr<TrialAdd> cloneShrPtr(Pair* pair,
                                   Criteria* criteria) const;

  // experimental double add if set to 1
  int addTwo = 0;

 protected:
  string molType_;   //!< type of molecule to add
  int molid_;        //!< index of molecule type

  /// attempt insertion
  void attempt1_();

  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(Pair *pair,
                                      Criteria *criteria) const;
};

/// Factory method
shared_ptr<TrialAdd> makeTrialAdd(Pair *pair, Criteria *criteria,
  const char* molType);

shared_ptr<TrialAdd> makeTrialAdd(const char* molType);

class MC;

/// Add a "TrialAdd" trial to the Monte Carlo object, mc for inserting moltype.
void addTrial(MC *mc, const char* moltype);

/// Add a "TrialAdd" trial to the Monte Carlo object, mc for inserting moltype.
void addTrial(shared_ptr<MC> mc, const char* moltype);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // TRIAL_ADD_H_

