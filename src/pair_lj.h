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
#include "./pair_1lrc.h"
#include "./functions.h"

namespace feasst {

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
class PairLJ : public PairLRC {
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
     *  cutType : method to handle the cut-off interactions (do not provide if molType none)
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

  /// Initialize Weeks-Chandler-Andersen for types itype and jtype.
  void initWCA(const int itype, const int jtype);

  // See overloaded virtual in base class.
  void initLRC();

  /// Initialize screened electrostatic interaction (Yukawa).
  /// \f$ U(r) = A e^{-K r}/r \f$
  void initScreenedElectro(const double A, const double K) {
    yukawa_ = 1; yukawaA_ = A; yukawaK_ = K;
  }

  // Similar as above, but sets yukawa flag
  void initScreenedElectro(const double A, const double K, const int yukawa) {
    yukawa_ = yukawa; yukawaA_ = A; yukawaK_ = K;
  }

  // set Aij for yukawa
  void initScreenedElectroIJ(const int itype, const int jtype, const double A,
                             const double K);

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

  vector<vector<double> > lambda() const { return lambda_; };       //!< DOI: 10.1021/ja802124e

 protected:
  double peLJ_;   //!< total potential energy from lennard-jones interactions
  double deLJ_;   //!< lennard jones potential energy change

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

}  // namespace feasst

#endif  // PAIR_LJ_H_
