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
  return feasst::str(coord_);
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

}  // namespace feasst
