
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

  /// Return the dot product of position vector with self.
  double dot_product(const Position &position) const;

  /// Return the squared distance of self from the origin.
  double squared_distance() const;

  /// Return the distance of self from the origin.
  double distance() const;

  /// Return coordinates as a string.
  std::string str() const;

  // HWH optimized only
  std::vector<double> * get_coord() { return &coord_; }

 private:
  std::vector<double> coord_;
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
