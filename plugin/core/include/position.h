
#ifndef FEASST_CORE_POSITION_H_
#define FEASST_CORE_POSITION_H_

#include <string>
#include <vector>

namespace feasst {

/**
  Positions are a set of coordinates (abbreviated here as "coord") in the
  Euclidean coordinate system. The number of coordinates is the dimensionality
  of the Euclidean space.
 */
class Position {
 public:
  /// Initialize a dimensional space.
  Position() {}

  /// Initialize coordinates by brace initialized position vector.
  explicit Position(std::vector<double> vec) { coord_ = vec; }

  /// Return a copy of the position vector.
  const std::vector<double>& coord() const { return coord_;}

  /// Set position vector.
  Position& set_vector(const std::vector<double> &doubleVec) {
    coord_ = doubleVec;
    return *this;
  }

  /// Get coordinate value of one dimension.
  double coord(const int dimension) const;

  /// Set coordinate value of one dimension.
  void set_coord(const int dimension, const double coord);

  /// Return the dimensionality of the position.
  int size() const { return coord_.size(); }

  /// Return the dimensionality of the position.
  int dimension() const { return size(); }

  /// Set the position of self to the origin in 3D space.
  /// HWH depreciate
  void set_to_origin_3D() { set_vector({0., 0., 0.}); }

  /// Set the position of self to the origin.
  void set_to_origin();

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
  bool is_equal(const Position& position) const;

  // HWH optimized only
  std::vector<double> * get_coord() { return &coord_; }

 private:
  std::vector<double> coord_;
};

/**
  Spherical coordinates are defined as follows:
  The first coordinate is rho >= 0. rho is the distance from origin.
  The second coordinate is theta.
  theta is the angle between x-axis and projection of line on x-y axis.
  The third and final coordinate is phi, 0 <= phi <= PI.
  phi is the angle between z-axis and line.
 */
class PositionSpherical : public Position {
 public:
  Position cartesian();
};

class SpatialEntity {
 public:
  /// Set the Position of the Site.
  void set_position(const Position& position) { position_ = position; }

  /// Add to the position.
  void add_position(const Position& position) { position_.add(position); }

  /// Return the Position of the Site.
  const Position& position() const { return position_; }

  /// Set the coordinate of one dimension of the Particle.
  void set_coordinate(const int dimension, const double coord) {
    position_.set_coord(dimension, coord);
  }

 private:
  Position position_;
};

}  // namespace feasst

#endif  // FEASST_CORE_POSITION_H_
