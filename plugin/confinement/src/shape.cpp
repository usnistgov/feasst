
#include <cmath>
#include <algorithm>
#include <memory>
#include "confinement/include/shape.h"
#include "utils/include/debug.h"

namespace feasst {

std::map<std::string, std::shared_ptr<Shape> >& Shape::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Shape> >* ans =
     new std::map<std::string, std::shared_ptr<Shape> >();
  return *ans;
}

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

void Shape::serialize(std::ostream& ostr) const { ERROR("not implemented"); }

std::shared_ptr<Shape> Shape::create(std::istream& istr) const {
  ERROR("not implemented");
}

std::shared_ptr<Shape> Shape::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr);
}

class MapShapeIntersect {
 public:
  MapShapeIntersect() {
    ShapeIntersect().deserialize_map()["ShapeIntersect"] = MakeShapeIntersect();
  }
};

static MapShapeIntersect mapper_ = MapShapeIntersect();

ShapeIntersect::ShapeIntersect(
    std::shared_ptr<Shape> shape1,
    std::shared_ptr<Shape> shape2) {
  shape1_ = shape1;
  shape2_ = shape2;
}

double ShapeIntersect::nearest_distance(const Position& point) const {
  const double dist1 = shape1_->nearest_distance(point),
               dist2 = shape2_->nearest_distance(point);
  if (dist1 > dist2) {
    return dist1;
  } else {
    return dist2;
  }
}

class MapShapeUnion {
 public:
  MapShapeUnion() {
    ShapeUnion().deserialize_map()["ShapeUnion"] = std::make_shared<ShapeUnion>();
  }
};

static MapShapeUnion mapper_shape_union_ = MapShapeUnion();

ShapeUnion::ShapeUnion(
    std::shared_ptr<Shape> shape1,
    std::shared_ptr<Shape> shape2) {
  shape1_ = shape1;
  shape2_ = shape2;
}

double ShapeUnion::nearest_distance(const Position& point) const {
  const double dist1 = shape1_->nearest_distance(point),
               dist2 = shape2_->nearest_distance(point);
  if (dist1 < dist2) {
    return dist1;
  } else {
    return dist2;
  }
}

}  // namespace feasst
