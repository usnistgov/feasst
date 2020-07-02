#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/site.h"
#include "configuration/include/model_params.h"
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

double ModelHardShape::energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) const {
  const int type = site.type();
  const double sigma = model_params.sigma().value(type);
  if (shape()->is_inside(wrapped_site, sigma)) {
    return 0.;
  }
  return NEAR_INFINITY;
}

}  // namespace feasst
