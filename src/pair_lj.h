/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_LJ_H_
#define PAIR_LJ_H_

#include <vector>
#include <map>
#include "./pair.h"
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
class PairLJ : public Pair {
 public:
  /// Constructor
  PairLJ(Space* space,
    /**
     * allowed string key pairs (e.g. dictionary):
     *
     *  molType : file which describes molecule
     *
     *  - (default): /path/to/feasst/forcefield/data.lj
     *
     *  - /path/to/lmp/data.file
     *
     *  - none: does not load a molecule, and therefore skips the rest of the
     *          present in this constructor, which now must be carried out
     *          manually. Providing other key pairs will result in an error.
     *          Warning: "none" without any further initialization of the Pair
     *          class may result in a segmentation fault.
     *
     *  molTypeInForcefield : file in Space.install_dir()/forcefield which
     *                        describes molecule. Only when molType is empty.
     *
     *  cutType : method to handle the cut-off interactions
     *
     *  - lrc(default): analytical long-range corrections
     *
     *  - cutShift: cut and shift energy to zero by a constant (no lrc)
     *
     *  - linearShift: cut and shift force to zero by a constant
     *                 and shift energy by a linear term (no lrc)
     *
     *  - none: do absolutely nothing about the cutoff (not recommended)
     */
    const argtype &args = argtype());

  /// Constructor
  PairLJ(shared_ptr<Space> space, const argtype &args = argtype())
    : PairLJ(space.get(), args) {}

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

  /// Initialize alpha parameter as a continuous parameter (non optimized).
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

  /**
   * Potential energy and forces of all particles.
   *  if flag == 0, dummy calculation
   *  if flag == 1, all config calculation
   */
  double allPartEnerForce(const int flag);

  void initForces(const int flag) {
    ASSERT(flag == 0 || expType_ == 0,
      "forces hard coded for expType_ == 0");
    Pair::initForces(flag);
  }

  // stores, restores or updates variables to avoid order recompute of entire
  //  configuration after every change
  void update(const vector<int> mpart, const int flag, const char* uptype);

  PairLJ(Space* space, const char* fileName);
  ~PairLJ() {}
  virtual PairLJ* clone(Space* space) const {
    PairLJ* p = new PairLJ(*this); p->reconstruct(space); return p;
  }

  /// Write restart file.
  virtual void writeRestart(const char* fileName);

  // read-only access of protected variables
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

  // See comments of derived class from Pair
  void pairSiteSite_(const int &iSiteType, const int &jSiteType, double * energy,
    double * force, int * neighbor, const double &dx, const double &dy,
    const double &dz);

  // Check for optimized loops
  virtual double pairLoopSite_(
    const vector<int> &siteList,
    const int noCell = 0);

  /// Overload default such that siteList is all sites in space
  double pairLoopSite_(const int noCell = 0) {
    const vector<int> &allSites = space_->listAtoms();
    return pairLoopSite_(allSites, noCell); }

  // defaults in constructor
  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairLJ> makePairLJ(Space* space,
  const argtype &args = argtype());

/// Factory method
shared_ptr<PairLJ> makePairLJ(shared_ptr<Space> space,
  const argtype &args = argtype());

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_LJ_H_
