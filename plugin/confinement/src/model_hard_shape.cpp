#include "confinement/include/model_hard_shape.h"

namespace feasst {

class MapModelHardShape {
 public:
  MapModelHardShape() {
    ModelHardShape().deserialize_map()["ModelHardShape"] =
      std::make_shared<ModelHardShape>();
  }
};

static MapModelHardShape map_model_hard_shape_ = MapModelHardShape();

void ModelHardShape::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  ShapedEntity::serialize(ostr);
  feasst_serialize_version(4276, ostr);
}

ModelHardShape::ModelHardShape(std::istream& istr)
  : ModelOneBody(), ShapedEntity(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 4276, "unrecognized verison: " << version);
}

}  // namespace feasst
