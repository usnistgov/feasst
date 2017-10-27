/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef TRIAL_CONFSWAP_OMP_H_
#define TRIAL_CONFSWAP_OMP_H_

#include <memory>
#include <string>
#include <vector>
#include "./trial.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Attempt to swap configurations with inter or intra processor stored state
 * as described in: http://dx.doi.org/10.1063/1.4918557 .
 */
class TrialConfSwapOMP : public Trial {
 public:
  /// Constructor
  TrialConfSwapOMP(Pair *pair, Criteria *criteria);

  /// Initialialize order parameter type.
  void initMType(const char* otype) { orderType_.assign(otype); }

  /// Add order parameter which overlaps with a given processor.
  void addProcOverlap(
    const double order,   //!< order parameter value
    TrialConfSwapOMP* trial,  //!< pointer to swap trial on other processor
    const double dbeta = 0.,   //!< change in beta of overlapping processor
    /// change in lnz of overlapping processor
    const double dlnz = 0.);

  /// Given order, return index (or -1 if not overlapping).
  int order2index(const double order);

  // functions for read-only access of private data-members
  vector<shared_ptr<Space> > confIntra() { return confIntra_; }
  vector<double> pe() { return pe_; }
  double orderTolerance() { return orderTolerance_; }
  vector<TrialConfSwapOMP*> trialSwapInter() { return trialSwapInter_; }

  /// This constructor is not often used, but its purpose is to initialize trial
  /// for interface before using reconstruct to set object pointers.
  TrialConfSwapOMP();

  /// Write restart file.
  void writeRestart(const char* fileName);

  /// Construct from restart file.
  TrialConfSwapOMP(const char* fileName, Pair *pair,
                   Criteria *criteria);
  ~TrialConfSwapOMP() {}
  TrialConfSwapOMP* clone(Pair* pair, Criteria* criteria) const {
    TrialConfSwapOMP* t = new TrialConfSwapOMP(*this);
    t->reconstruct(pair, criteria); return t; }
  shared_ptr<TrialConfSwapOMP> cloneShrPtr
    (Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialConfSwapOMP, Trial>
    (cloneImpl(pair, criteria))); }

 protected:
  string orderType_;    //!< obtain orderType from criteria upon initialization

  /// randomly attempt trial
  void attempt1_();

  /// list of order parameters that overlap with processors
  vector<double> order_;

  vector<double> pe_;   //!< potential energy of stored configuration

  /// given overlapping region, store configuration of current processor
  vector<shared_ptr<Space> > confIntra_;

  /// pointer to stored configurations of neighboring processor
  vector<TrialConfSwapOMP*> trialSwapInter_;

  vector<double> dbeta_;   //!< change in beta
  vector<double> dlnz_;    //!< change in lnz
  double orderTolerance_;  //!< tolerance in order2index

  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl
    (Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialConfSwapOMP> t = make_shared<TrialConfSwapOMP>(*this);
    t->reconstruct(pair, criteria); return t; }
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // TRIAL_CONFSWAP_OMP_H_

