/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef BASE_RANDOM_H_
#define BASE_RANDOM_H_

#include "./base.h"
#include "./random.h"
#include "./accumulator_vec.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * This class makes random number generation available to its derived classes
 * by storing a pointer to a random number generator, handing its
 * initialization and restarting, and providing a number of mathematical
 * operations for random number generation.
 */
class BaseRandom : public Base {
 public:
  BaseRandom();

  /// Construction from checkpoint file.
  explicit BaseRandom(const char* fileName);
  virtual ~BaseRandom() {};

  /// Write restart file for random number generator, if initialized.
  void writeRngRestart(const char* fileName);

  /// Initialize random number generator by pointer.
  void initRNG(shared_ptr<Random> ran);

  /** Initialize random number generator by seed. A zero seed signals use
   *  of the compiler random number generator to generate the seed.*/
  void initRNG(unsigned long long seed = 0);

  /// Initialize random number generator by checkpoint file.
  void initRNG(const char* fileName);

  /// Clear random number generator.
  void clearRNG();

  /// Return uiform random doubleprecision number between 0 and 1.
  double uniformRanNum();

  /// Return uniform random integer between range min and max, inclusive.
  int uniformRanNum(const int min, const int max);

  /// Return standard normal random number using Box-Muller method
  double stdNormRanNum();

  /// Multivariate Gaussian random number of given variance and dimension.
  vector<double> stdRanNum(const double variance, const int dimen);

  /// Return Gaussian distribution random number.
  double gaussRanNum(const double variance, const double av) {
    return av + variance*stdNormRanNum();
  }

  /** Return Gassuian distribution random number using implementation from
   *  Frenkel and Smit, "Understanding Molecular Simulation".
   *  Note that this implementation is slower than the Box-Muller method. */
  double gaussRanNumFS();

  /// Return a random hash of given length.
  std::string randomHash(const int length = 5);

  /** Return a random quaternion using the method of:
   *  Franz J. Vesely, J. Comput. Phys., 47, 291-296 (1982). */
  vector<double> quatRandom();

  /** Return a quaternion which is randomly perturbed from the unit vector
   *  (e.g., identity rotation matrix) by an amount "maxPerturb" */
  vector<double> quatRandom(const double maxPerturb);

  /** Return random position within spherical shell defined by rabove and rbelow
   * in dim dimensions */
  vector<double> ranShell(const double rabove, const double rbelow,
                          const int dim = 3);

  /// Return random position on unit sphere in dim dimensions.
  vector<double> ranUnitSphere(const int dim = 3);

  /** Given cumulative discrete probability distribution, cpdf,
   *  return chosen integer element from uniform probability distribution. */
  int ranFromCPDF(const vector<double> &cpdf);

  /** Return angle selected from probability distribution associated with
   *  harmonic bond bending energy, \f$ \beta U=k0(t-t0)^{power} \f$ from
   *  Frenkel and Smit, page 343, below Equation 13.3.6. */
  double ranAngle(const double k0, const double t0, const int power = 2);

  /** Return bond length selected from probability distribution associated with
   *  harmonic bond bending energy, \f$ \beta*U=k0*(l-l0)^{power} \f$. */
  double ranBond(const double k0, const double l0, const int power = 2);

  /**
   * Return random euler angles which yeild uniform sampling on a surface of a
   * sphere (e.g., uniform random orientation).
   *
   * Euler angles \f$(\phi, \theta, \psi)\f$ are defined by rotation
   * \f$\phi\f$ about z-axis, rotation \f$\theta\f$ about new x-xis,
   * then rotation \f$\psi\f$ about new z-axis.
   * The angles exist in the range \f$phi [-\pi,\pi], \theta [0, \pi],
   * \psi [-\pi, \pi]\f$.
   */
  vector<double> eulerRandom();

 protected:
  shared_ptr<Random> ranNum_;      //!< pointer to random number generator

  /// Derived objects may preform additional reconstruction.
  void reconstructDerived_();
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // BASE_RANDOM_H_

