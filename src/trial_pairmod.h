/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef TRIAL_PAIRMOD_H_
#define TRIAL_PAIRMOD_H_

#include <memory>
#include "./trial.h"

namespace feasst {

/**
 * Modify a parameter of the pair potential as the order parameter.
 */
class TrialPairMod : public Trial {
 public:
  /// Constructor
  TrialPairMod(Pair *pair, Criteria *criteria);

  /// This constructor is not often used, but its purpose is to initialize trial
  /// for interface before using reconstruct to set object pointers.
  TrialPairMod();

  /// ramp the chemical potential from value of min, located at
  // pair_->orderMin, to value of max for orderMax
  void initActivRamp(const double min, const double max) {
    lnzRampOn_ = 1; lnzMin_ = min; lnzMax_ = max;
  }

  // Overloaded from base class for status of specific trials.
  string printStat(const bool header = false);

  void writeRestart(const char* fileName);
  TrialPairMod(const char* fileName, Pair *pair,
               Criteria *criteria);
  ~TrialPairMod() {}
  TrialPairMod* clone(Pair* pair, Criteria* criteria) const {
    TrialPairMod* t = new TrialPairMod(*this);
    t->reconstruct(pair, criteria); return t;
  }
  shared_ptr<TrialPairMod> cloneShrPtr(
    Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialPairMod, Trial>(
      cloneImpl(pair, criteria)));
  }

 protected:
  /// store the ramped chemical potential values
  int lnzRampOn_;
  double lnzMin_, lnzMax_;

  void attempt1_();

  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(
    Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialPairMod> t = make_shared<TrialPairMod>(*this);
    t->reconstruct(pair, criteria);
    return t;
  }
};

/// Factory method
shared_ptr<TrialPairMod> makeTrialPairMod(Pair *pair,
  Criteria *criteria);

/// Factory method
shared_ptr<TrialPairMod> makeTrialPairMod();

class MC;
class WLTMMC;

void pairModTrial(MC *mc, const double maxMoveParam = -1);
void pairModTrial(shared_ptr<MC> mc, const double maxMoveParam = -1);
void pairModTrial(WLTMMC *mc, double maxMoveParam = -1);
void pairModTrial(shared_ptr<WLTMMC> mc, const double maxMoveParam = -1);

}  // namespace feasst

#endif  // TRIAL_PAIRMOD_H_

