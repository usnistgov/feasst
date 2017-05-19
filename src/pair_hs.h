/**
 * \file
 *
 * \brief hard sphere pair-wise interaction
 *
 */
#ifndef PAIR_HS_H_
#define PAIR_HS_H_

#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class PairHS : public Pair {
 public:
  PairHS(Space* space, const double rCut);
  PairHS(Space* space, const char* fileName);
  virtual ~PairHS() {}
  virtual PairHS* clone(Space* space) const {
    PairHS* p = new PairHS(*this); p->reconstruct(space); return p;
  }

  // defaults in constructor
  void defaultConstruction();

  /// write restart file
  virtual void writeRestart(const char* fileName);

  /// function to calculate forces, given positions
  virtual int initEnergy();

  /// potential energy of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag);
  void multiPartEnerAtomCutInner(const double &r2, const int &itype,
                                 const int &jtype);
  /// stores, restores or updates variables to avoid order recompute of entire
  //  configuration after every change
  void update(const vector<int> mpart, const int flag, const char* uptype);

 protected:
  double peSR_;
  double deSR_;
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_HS_H_

