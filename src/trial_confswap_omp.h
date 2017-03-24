/**
 * \file
 *
 * \brief
 *
 */

#ifndef TRIAL_CONFSWAP_OMP_H_
#define TRIAL_CONFSWAP_OMP_H_

#include "./trial.h"

class TrialConfSwapOMP : public Trial {
 public:
  TrialConfSwapOMP();
  TrialConfSwapOMP(Space *space, Pair *pair, Criteria *criteria);
  TrialConfSwapOMP(const char* fileName, Space *space, Pair *pair,
                   Criteria *criteria);
  ~TrialConfSwapOMP() {}
  TrialConfSwapOMP* clone(Space* space, Pair* pair, Criteria* criteria) const {
    TrialConfSwapOMP* t = new TrialConfSwapOMP(*this);
    t->reconstruct(space, pair, criteria); return t; }
  shared_ptr<TrialConfSwapOMP> cloneShrPtr
    (Space* space, Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialConfSwapOMP, Trial>
    (cloneImpl(space, pair, criteria))); }
  void writeRestart(const char* fileName);
  void defaultConstruction();

  /// randomly attempt trial
  void attempt1();

  /// add order parameter which overlaps with a given processor
  void addProcOverlap(const double order, TrialConfSwapOMP* trial)
    { addProcOverlap(order, trial, 0, 0); }
  void addProcOverlap(const double order, TrialConfSwapOMP* trial,
                      const double dbeta, const double dlnz);

  /// initialialize order parameter type
  void initMType(const char* otype) { orderType_.assign(otype); }

  /// given order, return index (or -1 if not overlapping)
  int order2index(const double order);

  // functions for read-only access of private data-members
  vector<shared_ptr<Space> > confIntra() { return confIntra_; }
  vector<double> pe() { return pe_; }
  double orderTolerance() { return orderTolerance_; }
  vector<TrialConfSwapOMP*> trialSwapInter() { return trialSwapInter_; }

 protected:
  string orderType_;    //!< obtain orderType from criteria upon initialization

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

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl
    (Space* space, Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialConfSwapOMP> t = make_shared<TrialConfSwapOMP>(*this);
    t->reconstruct(space, pair, criteria); return t; }
};

#endif  // TRIAL_CONFSWAP_OMP_H_

