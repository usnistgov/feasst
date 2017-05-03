/**
 * \file
 *
 * \brief pairwise interactions between hard circles in depletant
 */
#ifndef PAIR_HARD_CIRCLE_H_
#define PAIR_HARD_CIRCLE_H_

#include "./pair.h"

namespace feasst {

class PairHardCircle : public Pair {
 public:
  PairHardCircle(Space* space, const double rCut);
  PairHardCircle(Space* space, const char* fileName);
  virtual ~PairHardCircle() {}
  virtual PairHardCircle* clone(Space* space) const {
    PairHardCircle* p = new PairHardCircle(*this);
    p->reconstruct(space);
    return p;
  }

  // defaults in constructor
  void defaultConstruction();

  /// write restart file
  virtual void writeRestart(const char* fileName);

  virtual int initEnergy();   //!< function to calculate forces, given positions

  /// potential energy of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag);
  void multiPartEnerAtomCutInner(const double &r2, const int &itype,
    const int &jtype);

  /// potential energy and forces of all particles
  double allPartEnerForce(const int flag);

  /// potential energy and forces of all particles
  double allPartEnerForceNoCell();

  /// stores, restores or updates variables to avoid order recompute of entire
  //  configuration after every change
  void update(const vector<int> mpart, const int flag, const char* uptype);

  /// delete one particle
  void delPart(const int ipart) {delPartBase(ipart); }
  void delPart(const vector<int> mpart) {delPartBase(mpart); }

  //!< add one particle
  void addPart() {addPartBase(); }

  void initRDep(const double rDep) { rDep_ = rDep; }

  /// read-only access of protected variables
  double peSR() const { return peSR_; }
  double peSRone() const { return peSRone_; }

 protected:
  double peSR_;
  double deSR_;             //!< potential energy change
  double dCircle_;               //!< diameter of hard circle
  double rDep_;               //!< radius of depletant
};

}  // namespace feasst

#endif  // PAIR_HARD_CIRCLE_H_

