
#ifndef FEASST_MATH_POSITION_H_
#define FEASST_MATH_POSITION_H_

#include <string>
#include <vector>
#include "utils/include/arguments.h"
#include "math/include/constants.h"

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
  Position(const argtype& args);

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

  /**
    Set vector given Spherical coordinates (2 or 3-dimensional).
    Spherical coordinates are defined as follows:
    The first coordinate is rho >= 0. rho is the distance from origin.
    The second coordinate is theta.
    theta is the angle between x-axis and projection of line on x-y axis.
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
  double dot_product(const std::vector<double> vec) const;

  /// Return the cross product of position with self.
  Position cross_product(const Position& position) const;

  /// Return the squared distance of self from the origin.
  double squared_distance() const;

  /// Return the distance of self from the origin.
  double distance() const;

  /// Return coordinates as a string.
  std::string str() const;

  /// Return the cosine of the angle between self and given vector.
  double cosine(const Position& position) const;

  /// Normalize the position such that the distance from the origin is unity,
  /// but the direction from the origin is the same.
  void normalize();

  /// Return true if the given position is equal to self.
  bool is_equal(const Position& position,
                const double tolerance = NEAR_ZERO) const;

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
  Position(std::istream& istr);

 private:
  std::vector<double> coord_;
};

class SpatialEntity {
 public:
  SpatialEntity() {}

  /// Set the Position.
  void set_position(const Position& position) { position_ = position; }

  /// Add to the position.
  void add_position(const Position& position) { position_.add(position); }

  /// Return the Position.
  const Position& position() const { return position_; }

  /// Return the Position in a given dimension.
  const double position(const int dimension) const {
    return position_.coord(dimension); }

  /// Set the coordinate of one dimension.
  void set_coordinate(const int dimension, const double coord) {
    position_.set_coord(dimension, coord);
  }

  void serialize(std::ostream& ostr) const { position_.serialize(ostr); }
  SpatialEntity(std::istream& istr) { position_ = Position(istr); }

 private:
  Position position_;
};

}  // namespace feasst

#endif  // FEASST_MATH_POSITION_H_
