/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of this agreement (see LICENSE.txt) and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#ifndef PAIR_SQUAREWELL_H_
#define PAIR_SQUAREWELL_H_

#include <vector>
#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Square well interactions.
 * The hard distance is determined by sigma,
 * the square well distance is determined by rCut,
 * and the depth of the well is epsilon.
 */
class PairSquareWell : public Pair {
 public:
  /// Constructor
  PairSquareWell(Space* space, const double rCut);

  /// Initialize hard sphere interactions
  /// between particle types itype and jtype.
  void initHardSphere(const int itype, const int jtype);

  // function to calculate forces, given positions
  virtual void initEnergy();

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

  /// read-only access of protected variables
  double peSRone() const { return peSRone_; }

  PairSquareWell(Space* space, const char* fileName);
  virtual ~PairSquareWell() {}
  virtual PairSquareWell* clone(Space* space) const {
    PairSquareWell* p = new PairSquareWell(*this);
    p->reconstruct(space);
    return p;
  }

  /// write restart file
  virtual void writeRestart(const char* fileName);

 protected:
  double peSRone_;  //!< lennard jones potential energy from subset of particles
  double deSR_;     //!< lennard jones potential energy change

  // defaults in constructor
  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairSquareWell> makePairSquareWell(Space* space, const double rCut);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_SQUAREWELL_H_

