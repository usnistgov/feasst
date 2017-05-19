/**
 * \file
 *
 * \brief square well pair-wise interaction
 *
 * The hard distance is determined by sigma, and the square well distance is determined by rCut.
 * The depth of the well is epsilon
 */
#ifndef PAIR_SQUAREWELL_H_
#define PAIR_SQUAREWELL_H_

#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class PairSquareWell : public Pair {
 public:
  PairSquareWell(Space* space, const double rCut);
  PairSquareWell(Space* space, const char* fileName);
  virtual ~PairSquareWell() {}
  virtual PairSquareWell* clone(Space* space) const {
    PairSquareWell* p = new PairSquareWell(*this);
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

  /// potential energy and forces of all particles
  double allPartEnerForce(const int flag);
  double allPartEnerForceNoCell();

  /// potential energy and forces of all particles
  double allPartEnerForceCell();

  /// stores, restores or updates variables to avoid order recompute of entire
  //  configuration after every change
  void update(const vector<int> mpart, const int flag, const char* uptype);

  /// initialize hard sphere interactions between particle types itype and jtype
  void initHardSphere(const int itype, const int jtype);

  /// read-only access of protected variables
  double peSRone() const { return peSRone_; }

 protected:
  double peSRone_;  //!< lennard jones potential energy from subset of particles
  double deSR_;     //!< lennard jones potential energy change
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_SQUAREWELL_H_

