#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/site.h"
#include "configuration/include/model_params.h"
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

ModelLJShape::ModelLJShape(std::shared_ptr<Shape> shape,
  const argtype& args) : ModelOneBody(), ShapedEntity(shape) {
  args_.init(args);
  alpha_ = args_.key("alpha").dflt("3").dble();
}

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

double ModelLJShape::energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) const {
  const int type = site.type();
  const double sigma = model_params.sigma().value(type);
  const double epsilon = model_params.epsilon().value(type);
  const double distance = -shape()->nearest_distance(wrapped_site);
  if (distance <= 0) {
    return NEAR_INFINITY;
  } else {
    return epsilon * std::pow(distance/sigma, alpha_);
  }
}

}  // namespace feasst
