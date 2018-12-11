#include <cmath>
#include "core/include/position.h"
#include "core/include/debug.h"
#include "core/include/utils_io.h"

namespace feasst {

double Position::coord(const int dimension) const {
  ASSERT(dimension < size() && dimension >= 0, "size error");
  return coord_[dimension];
}

void Position::set_coord(const int dimension, const double coord) {
  ASSERT(dimension < size(), "size error");
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

}  // namespace feasst
