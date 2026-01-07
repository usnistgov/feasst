
#ifndef FEASST_MATH_CUBIC_SPLINE_H_
#define FEASST_MATH_CUBIC_SPLINE_H_

#include <map>
#include <string>
#include <vector>

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Linear spline is the same as linear interpolation.
  A piece-wise cubic polynomials with equal value at
  each knot point.

  \f$ y = C_0 + C_1 x \f$
 */
class LinearSpline1D {
 public:
  /// Initialize from a set y coordinates that are equally-spaced in x.
  /// The dimensionless variable z = (x - x_lower) / (x_upper - x_lower)
  /// spans the range from 0 to 1.
  explicit LinearSpline1D(const std::vector<double>& ycoord);

  /// Interpolate as a function of z.
  double interpolate(const double z) const;

  void serialize(std::ostream& ostr) const;
  explicit LinearSpline1D(std::istream& istr);
 private:
  std::vector<double> values_;
  double dz_;
};

typedef std::vector<std::vector<double> > vec2;

class LinearSpline2D {
 public:
  /// Initialize from a set y coordinates that are equally-spaced in x.
  /// The dimensionless variables xs = (x - x_lower) / (x_upper - x_lower)
  /// spans the range from 0 to 1.
  explicit LinearSpline2D(const vec2& ycoord);

  /// Interpolate as a function of scaled x and y.
  double interpolate(const double xs, const double ys) const;

  void serialize(std::ostream& ostr) const;
  explicit LinearSpline2D(std::istream& istr);
 private:
  vec2 values_;
  double dx_, dy_;
};

/**
  Quadratic spline is a piece-wise quadratic polynomials with equal value and
  equal derivative at each knot point.

  \f$ y = C_0 + C_1 x + C_2 x^2 \f$
 */
class QuadraticSpline1D {
 public:
  /// Initialize from a set y coordinates that are equally-spaced in x.
  /// The dimensionless variable z = (x - x_lower) / (x_upper - x_lower)
  /// spans the range from 0 to 1.
  explicit QuadraticSpline1D(const std::vector<double>& ycoord, const double right_endpoint_derivative=0.);

  /// Interpolate as a function of z.
  double interpolate(const double z) const;

  void serialize(std::ostream& ostr) const;
  explicit QuadraticSpline1D(std::istream& istr);
 private:
  std::vector<double> values_, derivs_;
  double dz_;
};

}  // namespace feasst

#endif  // FEASST_MATH_CUBIC_SPLINE_H_
