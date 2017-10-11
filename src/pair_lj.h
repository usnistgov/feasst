/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of this agreement (see LICENSE.txt) and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#ifndef PAIR_LJ_H_
#define PAIR_LJ_H_

#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Lennard Jones model, assuming one rCut and unit sigma and epsilon
 * \f$U_{LJ} = 4[r^{-12} - r^{-6}] \f$.
 */
class PairLJ : public Pair {
 public:
  /// Constructor
  PairLJ(Space* space,
    /// Interaction cut-off distance.
    const double rCut);

  virtual void initEnergy();   // function to calculate forces, given positions

  // potential energy of multiple particles
  virtual double multiPartEner(const vector<int> multiPart, const int flag);

  // potential energy of multiple particles optimized for no neighbor list
  //  updates
  double multiPartEnerNoNeigh(const vector<int> multiPart);

  // potential energy of multiple particles optimized for neighbor list updates
  double multiPartEnerNeigh(const vector<int> multiPart);

  // potential energy of multiple particles optimized for neighbor list updates
  virtual double multiPartEnerAtomCutNoNeigh(const vector<int> multiPart);

  // stores, restores or updates variables to avoid order recompute of entire
  //  configuration after every change
  void update(const vector<int> mpart, const int flag, const char* uptype);

  /**
   * Initialize cut and shifted potential.
   *  if flag == 0, do not shift
   *  if flag == 1, shift by potential value to zero at rCut
   */
  virtual void cutShift(const int flag = 0);

  /**
   * Initialize linear force shift potential.
   *  if flag == 0, do not shift
   *  if flag == 1, shift by potential and force value to zero at rCut
   */
  virtual void linearShift(const int flag = 0);

  // read-only access of protected variables
  double peLJ() const { return peLJ_; }
  double peLRC() const { return peLRC_; }
  double peLRCone() const { return peLRCone_; }

  /// Write restart file.
  virtual void writeRestart(const char* fileName);

  /// Construct with restart file.
  PairLJ(Space* space, const char* fileName);

  virtual ~PairLJ() {}
  virtual PairLJ* clone(Space* space) const {
    PairLJ* p = new PairLJ(*this); p->reconstruct(space); return p; }

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

/// Factor method
shared_ptr<PairLJ> makePairLJ(Space* space, const double rCut);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_LJ_H_

