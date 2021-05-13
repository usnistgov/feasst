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
  argtype args) : ModelLJShape(shape, &args) {
  check_all_used(args);
}

ModelLJShape::ModelLJShape(std::shared_ptr<Shape> shape,
  argtype * args) : ModelOneBody(), ShapedEntity(shape) {
  class_name_ = "ModelLJShape";
  alpha_ = dble("alpha", args, 3);
  disable_shift_ = boolean("disable_shift", args, false);
  shift_ = std::make_shared<ModelLJShapeEnergyAtCutoff>();
}

void ModelLJShape::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  ShapedEntity::serialize(ostr);
  feasst_serialize_version(1412, ostr);
  feasst_serialize(alpha_, ostr);
  feasst_serialize(disable_shift_, ostr);
  feasst_serialize_fstdr(shift_, ostr);
}

ModelLJShape::ModelLJShape(std::istream& istr)
  : ModelOneBody(istr), ShapedEntity(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1412, "unrecognized verison: " << version);
  feasst_deserialize(&alpha_, istr);
  feasst_deserialize(&disable_shift_, istr);
  // HWH for unknown reasons, this function template does not work
  // feasst_deserialize_fstdr(shift_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      shift_ = std::make_shared<ModelLJShapeEnergyAtCutoff>(istr);
    }
  }
}

double ModelLJShape::energy(const double epsilon,
              const double sigma,
              const double distance) const {
  TRACE("epsilon: " << epsilon);
  TRACE("sigma: " << sigma);
  TRACE("distance: " << distance);
  return epsilon * std::pow(sigma/distance, alpha_);
}

void ModelLJShape::precompute(const ModelParams& existing) {
  shift_->set_model(this); // note the model is used here for the computation
  shift_->set_param(existing);
  shift_->set_model(NULL); // remove model immediately
}

double ModelLJShape::energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) {
  const double distance = -shape()->nearest_distance(wrapped_site);
  TRACE("distance: " << distance);
  const int type = site.type();
  const double cutoff = model_params.cutoff().value(type);
  const double epsilon = model_params.epsilon().value(type);
  if (distance <= NEAR_ZERO && std::abs(epsilon) > NEAR_ZERO) {
    TRACE("here");
    TRACE(MAX_PRECISION << distance << " " << epsilon);
    return NEAR_INFINITY;
  } else if (distance >= cutoff) {
    return 0.;
  } else {
    const double sigma = model_params.sigma().value(type);
    const double en = energy(epsilon, sigma, distance);
    if (disable_shift_) {
      TRACE("en " << en);
      return en;
    } else {
      TRACE("shift " << shift_->value(type));
      TRACE("en " << en - shift_->value(type));
      return en - shift_->value(type);
    }
  }
}

}  // namespace feasst
