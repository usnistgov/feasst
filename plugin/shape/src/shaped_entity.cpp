#include "utils/include/serialize_extra.h"
#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "shape/include/shape.h"
#include "shape/include/shaped_entity.h"

namespace feasst {
  
ShapedEntity::ShapedEntity(std::shared_ptr<Shape> shape) { shape_ = shape; }

void ShapedEntity::set_shape(std::shared_ptr<Shape> shape) { shape_ = shape; }

const std::shared_ptr<Shape> ShapedEntity::shape() const { return shape_; }

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
