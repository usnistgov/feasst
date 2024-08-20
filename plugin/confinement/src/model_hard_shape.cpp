#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/constants.h"
#include "shape/include/shape_file.h"
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

ModelHardShape::ModelHardShape(argtype * args) : ModelHardShape() {
  set_shape(std::make_shared<ShapeFile>(args));
  cavity_ = boolean("cavity", args, true);
}
ModelHardShape::ModelHardShape(argtype args) : ModelHardShape(&args) {
  feasst_check_all_used(args);
}

void ModelHardShape::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  ShapedEntity::serialize(ostr);
  feasst_serialize_version(4276, ostr);
  feasst_serialize(cavity_, ostr);
}

ModelHardShape::ModelHardShape(std::istream& istr)
  : ModelOneBody(istr), ShapedEntity(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 4276, "unrecognized verison: " << version);
  feasst_deserialize(&cavity_, istr);
}

double ModelHardShape::energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) {
  TRACE("wrapped_site " << wrapped_site.str());
  const int type = site.type();
  TRACE("type " << type);
  const double sigma = model_params.select(sigma_index()).value(type);
  TRACE("sigma " << sigma);
  if (cavity_) {
    if (shape()->is_inside(wrapped_site, sigma)) {
      TRACE("inside cavity");
      return 0.;
    }
    TRACE("outside cavity");
  } else {
    if (!shape()->is_inside(wrapped_site, -sigma)) {
      TRACE("outside shape");
      return 0.;
    }
    TRACE("inside shape");
  }
  return NEAR_INFINITY;
}

}  // namespace feasst
