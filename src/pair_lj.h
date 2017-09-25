#ifndef PAIR_LJ_H_
#define PAIR_LJ_H_

#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * \file
 *
 * \brief lennard-jones pair-wise interaction
 *
 * The pair_lj class calculates pair-wise interactions with the Lennard-Jones classical potential function.
 * U is the potential energy, r is the separation distance between particles, eps is the energy-scale, sig is the length scale, and rCut is the cutoff of the interaction distance.
 *
 * U(r) = 4 eps ( (sig/r)^12 - (sig/r)^6 ) - U(rCut)
 * f_ = -grad(U)_ is the force
 * virial = sum( r_ . f_ ) is related to the pressure
 * where _ denotes vectors
 */
class PairLJ : public Pair {
 public:
  PairLJ(Space* space, const double rCut);
  PairLJ(Space* space, const char* fileName);
  virtual ~PairLJ() {}
  virtual PairLJ* clone(Space* space) const {
    PairLJ* p = new PairLJ(*this); p->reconstruct(space); return p; }

  /// write restart file
  virtual void writeRestart(const char* fileName);

  virtual void initEnergy();   //!< function to calculate forces, given positions

  /// potential energy of multiple particles
  virtual double multiPartEner(const vector<int> multiPart, const int flag);

  /// potential energy of multiple particles optimized for no neighbor list
  //  updates
  double multiPartEnerNoNeigh(const vector<int> multiPart);

  /// potential energy of multiple particles optimized for neighbor list updates
  double multiPartEnerNeigh(const vector<int> multiPart);

  /// potential energy of multiple particles optimized for neighbor list updates
  virtual double multiPartEnerAtomCutNoNeigh(const vector<int> multiPart);

  /// stores, restores or updates variables to avoid order recompute of entire
  //  configuration after every change
  void update(const vector<int> mpart, const int flag, const char* uptype);

  /// initialize cut and shifted potential
  virtual void cutShift(const int flag);

  /// initialize linear force shift potential
  virtual void linearShift(const int flag);

  /// read-only access of protected variables
  double peLJ() const { return peLJ_; }
  double peLRC() const { return peLRC_; }
  double peLRCone() const { return peLRCone_; }

 protected:
  double peLJ_;   //!< total potential energy from lennard-jones interactions
  double deLJ_;   //!< lennard jones potential energy change

  /// total potential energy from standard long range corrections
  double peLRC_;
  /// change in potential energy from standard long range corrections
  double deLRC_;

  /// potential energy from subset of particle standard long range corrections
  double peLRCone_;
  double peLRConePreCalc_;    //!< precalculated part of long range corrections
  double peShift_;          //!< shift potential energy by this amount
  bool cutShiftFlag_;       //!< flag to cut and shift potential by constant

  /// shift potential energy by this amount * (r-rc)
  double peLinearShift_;

  /// flag to cut and shift potential by linear term such that force=0 at rcut
  bool linearShiftFlag_;

  // defaults in constructor
  void defaultConstruction_();
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_LJ_H_

