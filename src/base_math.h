/**
 * \file
 *
 * \brief interface for all classes to inherit
 *
 */

#ifndef BASEMATH_H_
#define BASEMATH_H_

#include "./base_random.h"
#include "./accumulator_vec.h"

class BaseMath : public BaseRandom {
 public:
  BaseMath();
  virtual ~BaseMath() {}

  /// random quaterion
  vector<double> quatRandom();
  vector<double> quatRandom(const double maxDisp);

  /// random position within spherical shell defined by rabove and rebelow in
  //  dim dimensions
  vector<double> ranShell(const double rabove, const double rbelow,
                          const int dim);

  /// random position on unit sphere
  vector<double> ranUnitSphere(const int dim);

  /// given cumulative discrete probability distribution in vector of doubles,
  //  return chosen integer element from uniform probability distribution
  int ranFromCPDF(const vector<double> &cpdf);

  /// return angle selected from probability distribution associated with
  //  harmonic angle bending energy, \beta*U=k0*(t-t0)^pow
  double ranAngle(const double k0, const double t0, const int power = 2);

  /// return bond length selected from probability distribution associated
  //  with harmonic bond bending energy, \beta*U=k0*(l-l0)
  double ranBond(const double k0, const double l0, const int power = 2);

  /// return random euler angles which yeild uniform sampling on a surface of
  //  a sphere (e.g., uniform random orientation)
  vector<double> eulerRandom();

 protected:
};

#endif  // BASEMATH_H_

