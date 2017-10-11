/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of this agreement (see LICENSE.txt) and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#ifndef PAIR_LJ_MULTI_H_
#define PAIR_LJ_MULTI_H_

#include <vector>
#include "./pair_lj.h"
#include "./functions.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Lennard Jones particles with multiple sigmas, epsilons, rCuts, etc.
 * This class also includes other options, such as changing the exponential
 * parameter, alpha, \f$\alpha\f$.
 * There are also options to add Gaussians, Yukawa screened electrostatics,
 * Weeks-Chandler-Andersen and the LJ-lambda potential.
 * There are also long range corrections, cut and shift, and shifted-force cut
 * off methods.
 *
 * \f$U_{LJ} = 4\epsilon [ (\sigma/r)^{2\alpha} - (\sigma/r)^\alpha ] \f$.
 *
 */
class PairLJMulti : public PairLJ {
 public:
  /// Constructor
  /// @parameter rCut Interaction cut off distance.
  PairLJMulti(Space* space, const double rCut);

  /// Initialize cut and shifted potential for all particle types if flag == 1.
  /// \f$ U = U_{LJ} - U(rCut)\f$ when \f$r < rCut\f$.
  void cutShift(const int flag);

  /// Initialize linear force shift potential for all particle types when flag
  /// == 1.
  /// \f$ U = U_{LJ} - U(rCut) - (r-rCut)\left.\frac{\partial U}{\partial r}
  ///     \right|_{rCut}. \f$
  void linearShift(const int flag = 1);

  /// Initialize cut and shift for types itype and jtype if flag == 1.
  void cutShiftijset(const int itype, const int jtype, const int flag);

  /// Initialize linear force shift for types itype and jtype if flag == 1.
  void linearShiftijset(const int itype, const int jtype, const int flag);

  /// Initialize Weeks-Chandler-Andersen for types itype and jtype.
  void initWCA(const int itype, const int jtype);

  /// Initialize long-range corrections for all particle types.
  void initLRC();

  /// Initialize exponential type, alpha.
  //   type0 is 12-6: alpha=6
  //   type1 is 24-12: alpha=12
  //   type2 is 2alpha - alpha; alpha=16.6755
  //   type3 is 2alpha - alpha; alpha=50
  //   type4 is 2alpha - alpha; alpha=128
  //   type5 is 2alpha - alpha; alpha=24
  //   type6 is 2alpha - alpha; alpha=18
  void initExpType(const int type);

  /// Initialize alpha parameter as a continuous parameter (none optimized).
  void initAlpha(const double alpha = 12.);

  /// Initialize screened electrostatic interaction (Yukawa).
  /// \f$ U(r) = A e^{-K r}/r \f$
  void initScreenedElectro(const double A, const double K) {
    yukawa_ = 1; yukawaA_ = A; yukawaK_ = K;
  }

  // HWH: depreciated?
  void initScreenedElectro(const double A, const double K, const int yukawa) {
    yukawa_ = yukawa; yukawaA_ = A; yukawaK_ = K;
  }

  /// Set the order parameter value.
  void setOrder(const double order);

  /**
   * Add a gaussian on the potential.
   * \f$ U(r) = height * exp(-((r-position)/spread)^2) \f$
   */
  void addGaussian(const double height, const double position,
                   const double spread);

  /// Set lambda parameter for particle types iType and jType.
  /// http://dx.doi.org/10.1021/ja802124e
  void setLambdaij(const double iType, const double jType, const double lambda);

  /// Return the LRC contribution of one particle, ipart.
  double computeLRC(const int ipart);

  // read-only access to protected variables
  vector<double> rCutMax() const { return rCutMax_; }

  void initEnergy();     //!< function to calculate forces, given positions

  /// potential energy of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag);

  void multiPartEnerAtomCutInner(const double &r2, const int &itype,
                                 const int &jtype);

  /**
   * Potential energy and forces of all particles.
   *  if flag == 0, dummy calculation
   *  if flag == 1, all config calculation
   */
  double allPartEnerForce(const int flag);

  /// inner loop for potential energy and forces of all particles
  void allPartEnerForceInner(const double &r2, const double &dx,
    const double &dy, const double &dz, const int &itype, const int &jtype,
    const int &iMol, const int &jMol);

  PairLJMulti(Space* space, const char* fileName);
  ~PairLJMulti() {}
  virtual PairLJMulti* clone(Space* space) const {
    PairLJMulti* p = new PairLJMulti(*this); p->reconstruct(space); return p;
  }

  /// Write restart file.
  virtual void writeRestart(const char* fileName);

 protected:
  /// potential energy shift by constant for i-j type interactions
  vector<vector<double> > peShiftij_;

  /// potential energy shift by linear term for i-j type interactions
  vector<vector<double> > peLinearShiftij_;

  /// precalculation for tail corrections (long range corrections)
  vector<vector<double> > lrcPreCalc_;

  //!< exponential type. type0 is 12-6. type1 is 24-12. 2: 33.351-16.6755
  //  3: 100-50. 4: 256-128
  int expType_;

  int yukawa_;          //!< turn on yukawa interactions if 1
  double yukawaA_;      //!< U(r) = A*exp(-K*r)/r
  double yukawaK_;      //!< U(r) = A*exp(-K*r)/r
  double alpha_;        //!< exponential parameter

  // lambda parameters
  vector<vector<double> > lambda_;       //!< DOI: 10.1021/ja802124e
  int lambdaFlag_;      //!< default: 0 (off)

  // gaussian parameters
  int gaussian_;        //!< flag for guassian interacitons
  vector<vector<double> > gausParam_;

  // defaults in constructor
  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairLJMulti> makePairLJMulti(Space* space, const double rCut);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_LJ_MULTI_H_
