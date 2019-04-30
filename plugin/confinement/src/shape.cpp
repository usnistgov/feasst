
#include <cmath>
#include <algorithm>
#include <memory>
#include "confinement/include/shape.h"
#include "core/include/debug.h"

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
