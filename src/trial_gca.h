/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef TRIAL_GCA_H_
#define TRIAL_GCA_H_

#include <memory>
#include <vector>
#include "./trial.h"

namespace feasst {

/**
 * Rejection-Free Geometric Cluster Algorithm for Complex Fluids
 * Jiwen Liu and Erik Luijten
 * Phys. Rev. Lett. 92, 035504 – Published 23 January 2004
 */
class TrialGCA : public Trial {
 public:
  /// Constructor
  TrialGCA(Pair *pair, Criteria *criteria);

  /// This constructor is not often used, but its purpose is to initialize trial
  /// for interface before using reconstruct to set object pointers.
  TrialGCA();

  // tune parameters (e.g., based on acceptance)
  void tuneParameters();
  double targAcceptPer;      //!< target acceptance percentage

  /// Write restart file.
  void writeRestart(const char* fileName);

  /// Construct from restart file.
  TrialGCA(const char* fileName, Pair *pair, Criteria *criteria);
  ~TrialGCA() {}
  TrialGCA* clone(Pair* pair, Criteria* criteria) const {
    TrialGCA* t = new TrialGCA(*this); t->reconstruct(pair, criteria);
    return t; }
  shared_ptr<TrialGCA> cloneShrPtr(
    Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialGCA, Trial>(
      cloneImpl(pair, criteria)));
  }

 protected:
  void defaultConstruction_();

  void attempt1_();

  void attemptGCA_();

  /// recursively pivot molecule iMol
  void recursivePivot_(const vector<int> mpart, const vector<double> &xPivot,
                       vector<int> *pivotList, vector<int> *rejectNeigh,
                       vector<double> *rejectDe);

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(
    Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialGCA> t = make_shared<TrialGCA>(*this);
    t->reconstruct(pair, criteria);
    return t;
  }
};

/// Factory method
shared_ptr<TrialGCA> makeTrialGCA(Pair *pair, Criteria *criteria);

/// Factory method
shared_ptr<TrialGCA> makeTrialGCA();

class MC;
void gcaTrial(MC *mc, const int nMolTarg = -1);
void gcaTrial(shared_ptr<MC> mc, const int nMolTarg = -1);

}  // namespace feasst

#endif  // TRIAL_GCA_H_

