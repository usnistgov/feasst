/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of this agreement (see LICENSE.txt) and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#ifndef TRIAL_XSWAP_H_
#define TRIAL_XSWAP_H_

#include <memory>
#include "./trial.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Attempt to randomly swap the positions of two different particles.
 * To begin, randomly select two different particle types.
 * Then select a random particle of each type and attempt the swap.
 */
class TrialXSwap : public Trial {
 public:
  /// Constructor
  TrialXSwap(Pair *pair, Criteria *criteria);

  /// This constructor is not often used, but its purpose is to initialize trial
  /// for interface before using reconstruct to set object pointers.
  TrialXSwap();

  /// Construct from restart file.
  TrialXSwap(const char* fileName, Pair *pair,
             Criteria *criteria);
  ~TrialXSwap() {}
  TrialXSwap* clone(Pair* pair, Criteria* criteria) const {
    TrialXSwap* t = new TrialXSwap(*this);
    t->reconstruct(pair, criteria); return t;
  }
  shared_ptr<TrialXSwap> cloneShrPtr(
    Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialXSwap, Trial>(
      cloneImpl(pair, criteria)));
  }

 protected:
  void attempt1_();

  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(
    Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialXSwap> t = make_shared<TrialXSwap>(*this);
    t->reconstruct(pair, criteria);
    return t;
  }
};

/// Factory method
shared_ptr<TrialXSwap> makeTrialXSwap(Pair *pair,
  Criteria *criteria);

/// Factory method
shared_ptr<TrialXSwap> makeTrialXSwap();

class MC;
void xswapTrial(MC *mc);
void xswapTrial(shared_ptr<MC> mc);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // TRIAL_XSWAP_H_

