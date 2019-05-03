#include <cmath>
#include "core/include/position.h"
#include "core/include/debug.h"
#include "core/include/utils_io.h"
#include "core/include/constants.h"
#include "core/include/utils_math.h"

namespace feasst {

double Position::coord(const int dimension) const {
  return coord_[dimension];
}

void Position::set_coord(const int dimension, const double coord) {
  ASSERT(dimension < size(), "dimension(" << dimension << ") < size(" << size()
    << ")");
  coord_[dimension] = coord;
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
  return feasst_str(coord_);
}

void Position::set_to_origin() {
  for (double& value : coord_) {
    value = 0.;
  }
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

bool Position::is_equal(const Position& position) const {
  if (size() != position.size()) {
    return false;
  }
  for (int index = 0; index < static_cast<int>(size()); ++index) {
    if (std::abs(position.coord(index) - coord_[index]) > NEAR_ZERO) {
      return false;
    }
  }
  return true;
}

Position PositionSpherical::cartesian() {
  ASSERT(dimension() == 3, "requires 3D");
  const double rho = coord(0);
  const double theta = coord(1);
  const double phi = coord(2);
  const double sine_phi = sin(phi);
  Position pos;
  pos.set_vector({
    rho*sine_phi*cos(theta),
    rho*sine_phi*sin(theta),
    rho*cos(phi)});
  return pos;
}

void Position::multiply(const double constant) {
  for (double& coord : coord_) {
    coord *= constant;
  }
}

void Position::serialize(std::ostream& sstr) const {
  feasst_serialize_version(1, sstr);
  feasst_serialize(coord_, sstr);
}

Position::Position(std::istream& sstr) {
  feasst_deserialize_version(sstr);
  feasst_deserialize(&coord_, sstr);
}
  
Position Position::cross_product(const Position& position) const{
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

}  // namespace feasst
