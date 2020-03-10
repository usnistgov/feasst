#include <cmath>
#include "math/include/position.h"
#include "utils/include/debug.h"
#include "utils/include/utils_io.h"
#include "math/include/utils_math.h"

namespace feasst {

double Position::coord(const int dimension) const {
  return coord_[dimension];
}

void Position::set_coord(const int dimension, const double coord) {
  ASSERT(dimension < size(), "dimension(" << dimension << ") < size(" << size()
    << ")");
  coord_[dimension] = coord;
}

void Position::add_to_coord(const int dimension, const double coord) {
  ASSERT(dimension < size(), "dimension(" << dimension << ") < size(" << size()
    << ")");
  coord_[dimension] += coord;
}

void Position::add(const Position &position) {
  ASSERT(position.size() == size(), "size mismatch");
  for (int dim = 0; dim < size(); ++dim) {
    coord_[dim] += position.coord(dim);
  }
}

void Position::subtract(const Position &position) {
  ASSERT(position.size() == size(), "size mismatch");
  for (int dim = 0; dim < size(); ++dim) {
    coord_[dim] -= position.coord(dim);
  }
}

void Position::divide(const Position &position) {
  ASSERT(position.size() == size(), "size mismatch");
  for (int dim = 0; dim < size(); ++dim) {
    coord_[dim] /= position.coord(dim);
  }
}

double Position::dot_product(const Position &position) const {
  ASSERT(position.size() == size(), "size mismatch");
  double product = 0.;
  for (int dim = 0; dim < size(); ++dim) {
    product += position.coord(dim) * coord(dim);
  }
  return product;
}

double Position::dot_product(const std::vector<double> vec) const {
  Position pos;
  pos.set_vector(vec);
  return dot_product(pos);
}

double Position::squared_distance() const {
  double dist_sq = 0;
  for (const double& coord : coord_) {
    dist_sq += coord*coord;
  }
  return dist_sq;
}

double Position::distance() const {
  return std::sqrt(squared_distance());
}

std::string Position::str() const {
  return feasst_str(coord_, true);
}

void Position::set_to_origin() {
  for (double& value : coord_) {
    value = 0.;
  }
}

void Position::set_to_origin(const int dimension) {
  coord_.resize(dimension);
  set_to_origin();
}

double Position::cosine(const Position& position) const {
  return dot_product(position)/distance()/position.distance();
}

void Position::divide(const double denominator) {
  for (double& coord : coord_) {
    coord /= denominator;
  }
}

void Position::normalize() {
  divide(distance());
}

bool Position::is_equal(const Position& position,
                        const double tolerance) const {
  if (size() != position.size()) {
    return false;
  }
  for (int index = 0; index < static_cast<int>(size()); ++index) {
    if (std::abs(position.coord(index) - coord_[index]) > tolerance) {
      return false;
    }
  }
  return true;
}

Position& Position::set_from_spherical(const std::vector<double> &vec) {
  const int dimen = static_cast<int>(vec.size());
  if (dimen == 3) {
    set_from_spherical(vec[0], vec[1], vec[2]);
  } else if (dimen == 2) {
    set_from_spherical(vec[0], vec[1]);
  } else {
    ERROR("unrecognized dimension:" << vec.size());
  }
  return *this;
}

void Position::set_from_spherical(const double rho,
                                  const double theta,
                                  const double phi) {
  const double sine_phi = sin(phi);
  set_vector({
    rho*sine_phi*cos(theta),
    rho*sine_phi*sin(theta),
    rho*cos(phi)});
}

void Position::set_from_spherical(const double rho,
                                  const double theta) {
  set_vector({rho*cos(theta), rho*sin(theta)});
}

void Position::multiply(const double constant) {
  for (double& coord : coord_) {
    coord *= constant;
  }
}

void Position::serialize(std::ostream& sstr) const {
  feasst_serialize_version(3914, sstr);
  feasst_serialize(coord_, sstr);
}

Position::Position(std::istream& sstr) {
  const int version = feasst_deserialize_version(sstr);
  ASSERT(version == 3914, "unrecognized verison: " << version);
  feasst_deserialize(&coord_, sstr);
}

Position Position::cross_product(const Position& position) const {
  ASSERT(this->dimension() == 3 && position.dimension() == 3,
    "implemented for 3D only.");
  return Position().set_vector({
    this->coord(1) * position.coord(2) - this->coord(2) * position.coord(1),
    this->coord(2) * position.coord(0) - this->coord(0) * position.coord(2),
    this->coord(0) * position.coord(1) - this->coord(1) * position.coord(0)});
}

double Position::nearest_distance_to_axis(const Position& point1,
                                          const Position& point2) const {
  // point0 is self, pij = pi - pj
  Position p01 = *this, p02 = *this, p21 = point2;
  p01.subtract(point1);
  p02.subtract(point2);
  p21.subtract(point1);
  return ((p01.cross_product(p02)).distance())/p21.distance();
}

Position::Position(const argtype& args) {
  Arguments args_;
  args_.init(args);
  double x = 0, y = 0, z = 0;
  if (args_.key("x").used()) {
    x = args_.dble();
    if (args_.key("y").used()) {
      y = args_.dble();
      if (args_.key("z").used()) {
        z = args_.dble();
        set_vector({x, y, z});
      } else {
        set_vector({x, y});
      }
    } else {
      set_vector({x});
    }
  }
}

void Position::orthogonal(const Position& orthogonal) {
  ASSERT(orthogonal.size() == 3,
    "only implemented for 3D");
  if (size() != 3) set_to_origin(3);

  // find index with least absolute value
  int min = 0;
  if (std::abs(orthogonal.coord(1)) < std::abs(orthogonal.coord(min))) min = 1;
  if (std::abs(orthogonal.coord(2)) < std::abs(orthogonal.coord(min))) min = 2;

  // define orthogonal vector by setting the min to 0 and swapping the other
  // two indices, setting one of the swapped negative.
  if (min == 0) {
    coord_[0] = 0.;
    coord_[1] = -orthogonal.coord(2);
    coord_[2] = orthogonal.coord(1);
  } else if (min == 1) {
    coord_[0] = -orthogonal.coord(2);
    coord_[1] = 0.;
    coord_[2] = orthogonal.coord(0);
  } else {
    coord_[0] = -orthogonal.coord(1);
    coord_[1] = orthogonal.coord(0);
    coord_[2] = 0.;
  }
  normalize();
  ASSERT(std::abs(dot_product(orthogonal)) < 100*NEAR_ZERO,
    "self(" << str() << " is not orthogonal with: " << orthogonal.str() <<
    " if the dot product of the two is: " << dot_product(orthogonal));
}

Position& Position::set_vector(const std::vector<double> &vec) {
  coord_ = vec;
  return *this;
}

void Position::reflect(const Position& reflection_point) {
  for (int dim = 0; dim < dimension(); ++dim) {
    coord_[dim] = 2.*reflection_point.coord(dim) - coord_[dim];
  }
}

}  // namespace feasst
