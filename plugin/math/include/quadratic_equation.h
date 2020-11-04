
#ifndef FEASST_MATH_QUADRATIC_EQUATION_H_
#define FEASST_MATH_QUADRATIC_EQUATION_H_

#include <cmath>
#include "utils/include/debug.h"
#include "math/include/constants.h"

namespace feasst {

/// Solve the quadratic equation \f$ax^2+bx+c=0\f$
/// If the discriminant is <0, do not find the roots.
template<typename T>
void quadratic_equation(const T& a,
    const T& b,
    const T& c,
    T * discriminant,
    T * root1,
    T * root2) {
  *discriminant = b*b - 4*a*c;
  if (std::abs(a) < NEAR_ZERO) {
    *root1 = -c/b;
    return;
  }
  if (*discriminant < 0) return;
  const double sqrt_disc = std::sqrt(*discriminant);
  *root1 = 0.5*(-b + sqrt_disc)/a;
  *root2 = 0.5*(-b - sqrt_disc)/a;
}

}  // namespace feasst

#endif  // FEASST_MATH_QUADRATIC_EQUATION_H_
