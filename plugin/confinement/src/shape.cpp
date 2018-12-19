
#include <cmath>
#include <algorithm>
#include <memory>
#include "confinement/include/shape.h"
#include "core/include/debug.h"

namespace feasst {

bool Shape::is_inside(const Position& point) const {
  DEBUG(nearest_distance(point));
  if (nearest_distance(point) < 0) {
    return true;
  }
  return false;
}

bool Shape::is_inside(const Position& point, const double diameter) const {
  if (nearest_distance(point) + 0.5*diameter < 0) {
    return true;
  }
  return false;
}

HalfSpace& HalfSpace::set_direction(
    const double direction) {
  ASSERT(std::abs(direction) > 1e-15, "direction cannot be infinitesimal");
  if (direction > 0) {
    direction_ = 1.;
  } else {
    direction_ = -1.;
  }
  return *this;
}

double HalfSpace::nearest_distance(const Position& point) const {
  DEBUG("c " << point.coord(dimension_) << " " <<
                -1.*direction_*(point.coord(dimension_) - intersection_));
  return -1.*direction_*(point.coord(dimension_) - intersection_);
}

ShapeIntersect::ShapeIntersect(
    const std::shared_ptr<Shape> shape1,
    const std::shared_ptr<Shape> shape2) :
    shape1_(shape1),
    shape2_(shape2) {
}
double ShapeIntersect::nearest_distance(const Position& point) const {
  const double dist1 = shape1_->nearest_distance(point),
               dist2 = shape2_->nearest_distance(point);
  if (std::abs(dist1) < std::abs(dist2)) {
    return dist1;
  } else {
    return dist2;
  }
}

}  // namespace feasst
