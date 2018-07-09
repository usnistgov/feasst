/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_LRC_H_
#define PAIR_LRC_H_

#include <memory>
#include <vector>
#include "./pair.h"

namespace feasst {

/**
 * This class computes long-range corrections (LRC) due to the truncation of a
 * Lennard-Jones potential, as described on page 64 of Allen and Tildesley's
 * "Computer Simulation of Liquids" 1987.
 *
 * Special considerations must be made for the grand canonical ensemble,
 * where the number of particles changes. In this ensemble, it is difficult, if
 * not impossible, to add LRC contributions after the simulation is complete.
 * In addition, the contribution is nonlinear in number of particles.
 *
 * For multicomponent simulations,
 * \f$ U^{LRC}_{ij} = \sum_i \sum_j C_{ij} N_i N_j/V\f$,
 * where \f$\sum_i\f$ is a sum over particle types,
 * \f$N_i\f$ is the number of sites of type i, \f$r_c\f$ is the cut off,
 * \f$ C_{ij} = 2\pi\int_{r_c}^\infty U_{ij}g_{ij}r^2dr\f$, and
 * the radial distribution function, \f$g_{ij}\f$ is assumed to be unity.
 * \f$C_{ij}\f$ are initialized and assumed to remain constant during the
 * course fo the simulation (see below).
 *
 * When a site or multiple sites are added simultaneously,
 * \f$\Delta U^{LRC}_{ij} \sim N_i N_{ja} + N_{ia} N_j - N_{ia} N_{ja} \f$,
 * where \f$N_{ia}\f$ is the number of sites of type i that were added.
 */
class PairLRC : public Pair {
 public:
  /// Constructor
  PairLRC(Space* space, const argtype &args = argtype());

  /// Compute and store \f$C_{ij}\f$ for all particle types.
  virtual void initLRC();

  /// Return the long-range contribution.
  double computeLRC(
    /// compute contribution from subset of sites. If empty, consider all sites.
    const vector<int> msite = vector<int>());

  /// Initialize cut and shifted potential for all particle types if flag == 1.
  /// \f$ U = U_{LJ} - U(rCut)\f$ when \f$r < rCut\f$.
  virtual void cutShift(const int flag);

  /// Initialize cut and shift for types itype and jtype if flag == 1.
  void cutShiftijset(const int itype, const int jtype, const int flag);

  /// Initialize linear force shift for types itype and jtype if flag == 1.
  void linearShiftijset(const int itype, const int jtype, const int flag);

  /// Initialize linear force shift potential for all particle types when flag
  /// == 1.
  /// \f$ U = U_{LJ} - U(rCut) - (r-rCut)\left.\frac{\partial U}{\partial r}
  ///     \right|_{rCut}. \f$
  void linearShift(const int flag = 1);

  /// Initialize exponential type, alpha.
  //   type0 is 12-6: alpha=6
  //   type1 is 24-12: alpha=12
  //   type2 is 2alpha - alpha; alpha=16.6755
  //   type3 is 2alpha - alpha; alpha=50
  //   type4 is 2alpha - alpha; alpha=128
  //   type5 is 2alpha - alpha; alpha=24
  //   type6 is 2alpha - alpha; alpha=18
  void initExpType(const int type);

  /// Initialize alpha parameter as a continuous parameter (non optimized).
  void initAlpha(const double alpha = 12.);

  // read-only access of protected variables
  double peLRC() const { return peLRC_; }
  double peLRCone() const { return peLRCone_; }

  /// Write restart file.
  virtual void writeRestart(const char* fileName);

  PairLRC(Space* space, const char* fileName);

 protected:
  // defaults in constructor
  void defaultConstruction_();

  /// total potential energy from standard long range corrections
  double peLRC_;
  /// change in potential energy from standard long range corrections
  double deLRC_;

  /// potential energy from subset of particle standard long range corrections
  double peLRCone_;
  double peLRConePreCalc_;    //!< precalculated part of long range corrections

  /// precalculation for tail corrections (long range corrections)
  vector<vector<double> > lrcPreCalc_;

  double peShift_;          //!< shift potential energy by this amount
  bool cutShiftFlag_;       //!< flag to cut and shift potential by constant

  /// shift potential energy by this amount * (r-rc)
  double peLinearShift_;

  /// flag to cut and shift potential by linear term such that force=0 at rcut
  bool linearShiftFlag_;

  /// potential energy shift by constant for i-j type interactions
  vector<vector<double> > peShiftij_;

  //!< exponential type. type0 is 12-6. type1 is 24-12. 2: 33.351-16.6755
  //  3: 100-50. 4: 256-128
  int expType_;

  int yukawa_;          //!< turn on yukawa interactions if 1
  double yukawaA_;      //!< U(r) = A*exp(-K*r)/r
  vector<vector<double> > yukawaAij_;      //!< U(r) = A*exp(-K*r)/r
  double yukawaK_;      //!< U(r) = A*exp(-K*r)/r
  vector<vector<double> > yukawaKij_;      //!< U(r) = A*exp(-K*r)/r
  double alpha_;        //!< exponential parameter

  /// potential energy shift by linear term for i-j type interactions
  vector<vector<double> > peLinearShiftij_;
};

}  // namespace feasst

#endif  // PAIR_LRC_H_
