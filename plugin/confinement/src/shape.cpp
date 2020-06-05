#include "utils/include/serialize.h"
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

void Shape::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<Shape> Shape::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<Shape> Shape::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr);
}

void ShapedEntity::serialize(std::ostream& ostr) const {
  feasst_serialize_version(9249, ostr);
  feasst_serialize_fstdr(shape_, ostr);
}

ShapedEntity::ShapedEntity(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9249, "unrecognized verison: " << version);
  // feasst_deserialize_fstdr(shape_, istr);
  { // HWH for unknown reasons the above template function does not work
    int existing;
    istr >> existing;
    if (existing != 0) {
      shape_ = shape_->deserialize(istr);
    }
  }
}

}  // namespace feasst
