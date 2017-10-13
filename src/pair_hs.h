/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#ifndef PAIR_HS_H_
#define PAIR_HS_H_

#include <vector>
#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Hard sphere pair-wise interaction.
 * The diameter is determined by the sigma parameter described in Pair.
 */
class PairHS : public Pair {
 public:
  /// Constructor
  /// @param rCut interaciton cut off distance
  PairHS(Space* space, const double rCut);

  /// Initialize interactions.
  virtual void initEnergy();

  /// potential energy of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag);
  void multiPartEnerAtomCutInner(const double &r2, const int &itype,
                                 const int &jtype);
  /// stores, restores or updates variables to avoid order recompute of entire
  //  configuration after every change
  void update(const vector<int> mpart, const int flag, const char* uptype);

  /// potential energy and forces of all particles
  double allPartEnerForce(const int flag);

  /// inner loop for potential energy and forces of all particles
  void allPartEnerForceInner(const double &r2, const double &dx,
    const double &dy, const double &dz, const int &itype, const int &jtype,
    const int &iMol, const int &jMol);

  PairHS(Space* space, const char* fileName);
  virtual ~PairHS() {}
  virtual PairHS* clone(Space* space) const {
    PairHS* p = new PairHS(*this); p->reconstruct(space); return p;
  }

  /// write restart file
  virtual void writeRestart(const char* fileName);

 protected:
  double peSR_;
  double deSR_;

  // defaults in constructor
  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairHS> makePairHS(Space* space, const double rCut);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_HS_H_

