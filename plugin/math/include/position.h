
#ifndef FEASST_MATH_POSITION_H_
#define FEASST_MATH_POSITION_H_

#include <string>
#include <vector>
#include <memory>
#include "utils/include/arguments.h"

namespace feasst {

/**
  Positions are a set of coordinates (abbreviated here as "coord") in the
  Euclidean/Cartesian coordinate system.
  The number of coordinates is the dimensionality of the Euclidean space.
 */
class Position {
 public:
  Position() {}

  /**
    args:
    - x: x-coordinate
    - y: y-coordinate. Requires explicit x.
    - z: z-coordinate. Requires explicit y.
   */
  explicit Position(argtype args);

  /// Initialize coordinates by brace initialized position vector.
  explicit Position(std::vector<double> vec) { coord_ = vec; }

  /// Initialize coordinates on origin with given dimensionality.
  explicit Position(const int dimension) { set_to_origin(dimension); }

  /// Return a copy of the position vector.
  const std::vector<double>& coord() const { return coord_;}

  /// Set position vector.
  Position& set_vector(const std::vector<double> &vec);

  /// Set vector given Cartesian coordinates (same as above).
  Position& set_from_cartesian(const std::vector<double> &vec) {
    return set_vector(vec); }

  /// Add another coordinate dimension.
  void push_back(const double coord) { coord_.push_back(coord); }

  /**
    Set vector given Spherical coordinates (2 or 3-dimensional).
    See https://mathworld.wolfram.com/SphericalCoordinates.html
    Spherical coordinates are defined as follows:
    The first coordinate is rho >= 0. rho is the distance from origin.
    The second coordinate is theta.
    theta is the angle between x-axis and projection of the vector on x-y plane.
    In 3D, the third and final coordinate is phi, 0 <= phi <= PI.
    phi is the angle between z-axis and line.
   */
  Position& set_from_spherical(const std::vector<double> &vec);

  /// Same as above for 3D.
  void set_from_spherical(const double rho,
                          const double theta,
                          const double phi);

  /// Same as above for 2D.
  void set_from_spherical(const double rho, const double theta);

  /// Return the spherical coordinates (r, theta, phi).
  Position spherical() const;

  /// Get coordinate value of one dimension.
  double coord(const int dimension) const;

  /// Set coordinate value of one dimension.
  void set_coord(const int dimension, const double coord);

  /// Add to coordinate value of one dimension.
  void add_to_coord(const int dimension, const double coord);

  /// Return the dimensionality of the position.
  int size() const { return coord_.size(); }

  /// Return the dimensionality of the position.
  int dimension() const { return size(); }

  /// Set the position of self to the origin in 3D space.
  /// HWH depreciate
  void set_to_origin_3D() { set_vector({0., 0., 0.}); }

  /// Set the position of self to the origin.
  void set_to_origin();

  /// Resize to given dimension, then set to the origin.
  void set_to_origin(const int dimension);

  /// Add position vector to self.
  void add(const Position &position);

  /// Subtract the position vector from self.
  void subtract(const Position &position);

  /// Divide self by the position vector.
  void divide(const Position &position);

  /// Divide self by a constant.
  void divide(const double denominator);

  /// Multiply self by a constant.
  void multiply(const double constant);

  /// Return the dot product of position vector with self.
  double dot_product(const Position &position) const;

  /// Same as above, but with a vector.
  double dot_product(const std::vector<double>& vec) const;

  /// Return the cross product of position with self.
  Position cross_product(const Position& position) const;

  /// Return the squared distance of self from the origin.
  double squared_distance() const;

  /// Return the distance of self from the origin.
  double distance() const;

  /// Return the squared distance between self and position.
  double squared_distance(const Position& position) const;

  /// Return distance between self and position.
  double distance(const Position& position) const;

  /// Return coordinates as a string.
  std::string str() const;

  /// Return the cosine of the angle between self and given vector.
  double cosine(const Position& position) const;

  /// Normalize the position such that the distance from the origin is unity,
  /// but the direction from the origin is the same.
  void normalize();

  /// Return true if the given position is equal to self within tolerance.
  bool is_equal(const Position& position,
                const double tolerance) const;

  /// Same as above, but with a default tolerance of NEAR_ZERO.
  bool is_equal(const Position& position) const;

  /// Nearest distance to axis defined by two points.
  /// see http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  double nearest_distance_to_axis(const Position& point1,
                                  const Position& point2) const;

  /// Set self orthogonal to given position.
  void orthogonal(const Position& orthogonal);

  /// Reflect self about a point.
  void reflect(const Position& reflection_point);

  // HWH optimized only
  std::vector<double> * get_coord() { return &coord_; }

  void serialize(std::ostream& ostr) const;
  explicit Position(std::istream& istr);

 private:
  std::vector<double> coord_;
};

inline std::shared_ptr<Position> MakePosition(argtype args = argtype()) {
  return std::make_shared<Position>(args);
}

inline std::shared_ptr<Position> MakePosition(std::vector<double> vec) {
  return std::make_shared<Position>(vec);
}

inline std::shared_ptr<Position> MakePosition(const int dimension) {
  return std::make_shared<Position>(dimension);
}

}  // namespace feasst

#endif  // FEASST_MATH_POSITION_H_
