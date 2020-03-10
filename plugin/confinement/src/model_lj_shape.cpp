#include "confinement/include/model_lj_shape.h"

namespace feasst {

class MapModelLJShape {
 public:
  MapModelLJShape() {
    ModelLJShape().deserialize_map()["ModelLJShape"] =
      std::make_shared<ModelLJShape>();
  }
};

static MapModelLJShape map_model_hard_shape_ = MapModelLJShape();

void ModelLJShape::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  ShapedEntity::serialize(ostr);
  feasst_serialize_version(1412, ostr);
  feasst_serialize(alpha_, ostr);
}

ModelLJShape::ModelLJShape(std::istream& istr)
  : ModelOneBody(), ShapedEntity(istr) {
  const int version = feasst_deserialize_version(istr);
  feasst_deserialize(&alpha_, istr);
  ASSERT(version == 1412, "unrecognized verison: " << version);
}

}  // namespace feasst
