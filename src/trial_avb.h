/**
 * \file
 *
 * \brief trial moves for Monte Carlo
 *
 */

#ifndef TRIAL_AVB_H_
#define TRIAL_AVB_H_

#include "./trial.h"

class TrialAVB : public Trial {
 public:
  TrialAVB(const double pBias, const double rAbove,
           const double rBelow, const int avbType);
  TrialAVB(Space *space, Pair *pair, Criteria *criteria, const double pBias,
           const double rAbove, const double rBelow, const int avbType);
  TrialAVB(const char* fileName, Space *space, Pair *pair, Criteria *criteria);
  ~TrialAVB() {}
  TrialAVB* clone(Space* space, Pair *pair, Criteria *criteria) const
    { TrialAVB* t = new TrialAVB(*this); t->reconstruct(space, pair, criteria);
      return t; }
  shared_ptr<TrialAVB> cloneShrPtr(Space* space, Pair* pair,
    Criteria* criteria) const {
      return(std::static_pointer_cast<TrialAVB, Trial>
      (cloneImpl(space, pair, criteria))); }

  void attempt1() {}

  /// attempt aggregation-volume-bias move
  ///  Chen and Siepmann JPC B Vol 104 No 36 (2000)
  ///                    JPC B Vol 105 p11275 (2001)
  void avb1() {}
  void avb2() {}
  void avb3() {}

 protected:
  double pBias_;  //!< bias probability
  int avbType_;   //!< type of avb move

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl
    (Space* space, Pair *pair, Criteria *criteria) const {
      shared_ptr<TrialAVB> t = make_shared<TrialAVB>(*this);
      t->reconstruct(space, pair, criteria); return t; }
};

#endif  // TRIAL_AVB_H_

