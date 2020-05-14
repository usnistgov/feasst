/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef TRIAL_BETA_H_
#define TRIAL_BETA_H_

#include <memory>
#include "./trial.h"

namespace feasst {

/**
 * Attempt to change beta, the inverse temperature.
 *
 * There are a number of parameters to optimize. They depend on the system
 * but we don't have a very rigorous way to determine these. I would determine
 * them in the order below.
 *
 * 1. macrostate range (e.g., temperature)
 *
 * 2. bin width
 *
 * 3. MC trial weight to move between bins (e.g., temperature trial move)
 *
 * For a) I would shorten the range as much as possible while capturing the
 * transition of interest. Check if the shorter range simulation compares
 * reasonably with the larger-range, slower-converging simulation. If the
 * simulation does not sample beyond the transition of interest and it is
 * metastable, you could see the transition diminish and the two simulations
 * will not match (in which case I would more trust the larger-range).
 *
 * For b) it depends what resolution you want your data in the macrostate and
 * also what is the acceptance probability to transition. Probably there is an
 * optimal bin spacing similar to replica exchange MD temperature choice based
 * on overlapping distributions.
 *
 * For c) and based on your choice of a&b) if the MC trial is an inexpensive
 * move (temperature change is "free" based on system energy), you can do this
 * move frequently. But if you do it too frequently the simulation artificially
 * converges before sampling different configurations. It is safer to not
 * perform these moves too frequently, but infrequent moves slows the
 * convergence of the simulations.
 */
class TrialBeta : public Trial {
 public:
  /// Constructor
  TrialBeta(Pair *pair, Criteria *criteria);

  /// This constructor is not often used, but its purpose is to initialize trial
  /// for interface before using reconstruct to set object pointers.
  TrialBeta();

  // Overloaded from base class for status of specific trials.
  string printStat(const bool header = false);

  // Construct from restart file.
  TrialBeta(const char* fileName, Pair *pair, Criteria *criteria);
  ~TrialBeta() {}
  TrialBeta* clone(Pair* pair, Criteria* criteria) const {
    TrialBeta* t = new TrialBeta(*this); t->reconstruct(pair, criteria);
    return t; }
  shared_ptr<TrialBeta> cloneShrPtr(Pair* pair,
                                    Criteria* criteria) const {
    return(std::static_pointer_cast<TrialBeta, Trial>
      (cloneImpl(pair, criteria))); }

 protected:
  /// Implementation of trial.
  void attempt1_();

  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl
    (Pair *pair, Criteria *criteria) const {
      shared_ptr<TrialBeta> t = make_shared<TrialBeta>(*this);
      t->reconstruct(pair, criteria);
      return t;}
};

shared_ptr<TrialBeta> makeTrialBeta(Pair *pair,
  Criteria *criteria);

shared_ptr<TrialBeta> makeTrialBeta();

class MC;
class WLTMMC;

/// Add a "TrialBeta" object to the Monte Carlo object, mc
/// @param maxMoveParam optionally set max move, otherwise default value if -1.
void betaTrial(MC *mc, const double maxMoveParam = -1);

/// Add a "TrialBeta" object to the Monte Carlo object, mc
/// @param maxMoveParam optionally set max move, otherwise default value if -1.
void betaTrial(shared_ptr<MC> mc, const double maxMoveParam = -1);

/// Add a "TrialBeta" object to the flat histogram Monte Carlo object, mc
/// @param maxMoveParam optionally set max move.
///   This time, the default maxMoveParameter is set equal to the TMWL bin size
void betaTrial(WLTMMC *mc, double maxMoveParam = -1);

/// Add a "TrialBeta" object to the flat histogram Monte Carlo object, mc
/// @param maxMoveParam optionally set max move.
///   This time, the default maxMoveParameter is set equal to the TMWL bin size
void betaTrial(shared_ptr<WLTMMC> mc, const double maxMoveParam = -1);

}  // namespace feasst

#endif  // TRIAL_BETA_H_

