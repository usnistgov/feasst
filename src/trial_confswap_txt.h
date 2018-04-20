/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef TRIAL_CONFSWAPTXT_H_
#define TRIAL_CONFSWAPTXT_H_

#include "./trial.h"

namespace feasst {

/**
 *  HWH: Note: This class is depreciated but remains as a stub file for
 *  HWH: Note:  eventual support.
 *  Attempts configuration swaps between processors.
 *  These swaps may occur with variable beta and lnz, but not variable order
 *  parameter (e.g. nMol or pairOrder)
 */
class TrialConfSwapTXT : public Trial {
 public:
  TrialConfSwapTXT();
  TrialConfSwapTXT(Pair *pair, Criteria *criteria);
  TrialConfSwapTXT(const char* fileName, Pair *pair,
                   Criteria *criteria);
  ~TrialConfSwapTXT() {}
  TrialConfSwapTXT* clone(Pair* pair, Criteria* criteria) const {
    TrialConfSwapTXT* t = new TrialConfSwapTXT(*this);
    t->reconstruct(pair, criteria); return t; }
  shared_ptr<TrialConfSwapTXT> cloneShrPtr
    (Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialConfSwapTXT, Trial>
    (cloneImpl(pair, criteria))); }
  void writeRestart(const char* fileName);
  void defaultConstruction_();

  // add order parameter which overlaps with a given processor
  void addProcOverlap(const double order, const int proc) {
    addProcOverlap(order, proc, 0, 0); }
  void addProcOverlap(const double order, const int proc, const double dbeta,
                      const double dlnz);

  // initialize process number
  void initProc(const int proc) { proc_ = proc; }

  // initialialize order parameter type
  void initMType(const char* otype) { orderType_.assign(otype); }

  // functions for read-only access of private data-members
  int nOverlaps() const { return order_.size(); }
  int proc() const { return proc_; }

 protected:
  string orderType_;   //!< obtain orderType from criteria upon initialization
  int proc_;           //!< process number of current trial

  /**
   * randomly to swap configurations with inter or intra processor stored state
   */
  void attempt1_();

  /// order parameter in overlapping reigion (e.g. nMol or pairOrder)
  vector<double> order_;

  vector<int> process_;   //!< process number of neighboring processor
  vector<int> nLines_;    //!< number of lines in restart file
  vector<double> dbeta_;  //!< change in beta
  vector<double> dlnz_;   //!< change in lnz

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(
    Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialConfSwapTXT> t = make_shared<TrialConfSwapTXT>(*this);
    t->reconstruct(pair, criteria); return t;
  }
};

}  // namespace feasst

#endif  // TRIAL_CONFSWAPTXT_H_

