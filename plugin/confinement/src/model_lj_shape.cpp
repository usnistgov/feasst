#include <cmath>
#include "utils/include/max_precision.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "shape/include/shape_file.h"
#include "configuration/include/site.h"
#include "configuration/include/model_params.h"
#include "confinement/include/model_lj_shape.h"

namespace feasst {

FEASST_MAPPER(ModelLJShape,);

void ModelLJShape::parse_args_(argtype * args) {
  class_name_ = "ModelLJShape";
  alpha_ = dble("alpha", args, 3.);
  delta_ = dble("delta", args, 0.);
  disable_shift_ = boolean("disable_shift", args, false);
  shift_ = std::make_shared<ModelLJShapeEnergyAtCutoff>();
  wall_sigma_ = dble("wall_sigma", args, 0.);
  wall_epsilon_ = dble("wall_epsilon", args, 0.);
}

ModelLJShape::ModelLJShape(argtype * args) : ModelLJShape() {
  set_shape(std::make_shared<ShapeFile>(args));
  parse_args_(args);
}
ModelLJShape::ModelLJShape(argtype args) : ModelLJShape(&args) {
  feasst_check_all_used(args);
}

ModelLJShape::ModelLJShape(std::shared_ptr<Shape> shape,
  argtype args) : ModelLJShape(shape, &args) {
  feasst_check_all_used(args);
}

ModelLJShape::ModelLJShape(std::shared_ptr<Shape> shape,
  argtype * args) : ModelOneBody(), ShapedEntity(shape) {
  parse_args_(args);
}

void ModelLJShape::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  ShapedEntity::serialize(ostr);
  feasst_serialize_version(1412, ostr);
  feasst_serialize(alpha_, ostr);
  feasst_serialize(delta_, ostr);
  feasst_serialize(disable_shift_, ostr);
  feasst_serialize_fstdr(shift_, ostr);
  feasst_serialize(wall_sigma_, ostr);
  feasst_serialize_fstobj(mixed_sigma_, ostr);
  feasst_serialize(wall_epsilon_, ostr);
  feasst_serialize_fstobj(mixed_epsilon_, ostr);
}

ModelLJShape::ModelLJShape(std::istream& istr)
  : ModelOneBody(istr), ShapedEntity(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1412, "unrecognized verison: " << version);
  feasst_deserialize(&alpha_, istr);
  feasst_deserialize(&delta_, istr);
  feasst_deserialize(&disable_shift_, istr);
  // HWH for unknown reasons, this function template does not work
  // feasst_deserialize_fstdr(shift_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      shift_ = std::make_shared<ModelLJShapeEnergyAtCutoff>(istr);
    }
  }
  feasst_deserialize(&wall_sigma_, istr);
  feasst_deserialize_fstobj(&mixed_sigma_, istr);
  feasst_deserialize(&wall_epsilon_, istr);
  feasst_deserialize_fstobj(&mixed_epsilon_, istr);
}

double ModelLJShape::energy(const double epsilon,
              const double sigma,
              const double distance) const {
  TRACE("epsilon: " << epsilon);
  TRACE("sigma: " << sigma);
  TRACE("distance: " << distance);
  return epsilon * std::pow(sigma/(distance + delta_), alpha_);
}

void ModelLJShape::precompute(const ModelParams& existing) {
  ModelOneBody::precompute(existing);
  if (std::abs(wall_sigma_) > NEAR_ZERO && mixed_sigma_.size() == 0) {
    //mixed_sigma_ = Sigma();  // reset in case multiple precompute
    const ModelParam& fluid_sig = existing.select("sigma");
    for (int type = 0; type < existing.size(); ++type) {
      mixed_sigma_.add(0.5*(fluid_sig.value(type) + wall_sigma_));
    }
    ASSERT(mixed_sigma_.size() == fluid_sig.size(), "error");
  }
  if (std::abs(wall_epsilon_) > NEAR_ZERO && mixed_epsilon_.size() == 0) {
    //mixed_epsilon_ = Epsilon();  // reset in case multiple precompute
    const ModelParam& fluid_eps = existing.select("epsilon");
    for (int type = 0; type < existing.size(); ++type) {
      const double walle = wall_epsilon_;
      const double fluide = fluid_eps.value(type);
      DEBUG("fluid eps " << fluide);
      DEBUG("wall_epsilon_ " << walle);
      mixed_epsilon_.add(feasst::sgn(walle)*std::sqrt(fluide*std::abs(walle)));
    }
    DEBUG("mixed eps " << mixed_epsilon_.str());
    ASSERT(mixed_epsilon_.size() == fluid_eps.size(), "error");
  }
  // compute shift after possible mixing with all parameters
  shift_->set_model(this); // note the model is used here for the computation
  shift_->set_param(existing);
  shift_->set_model(NULL); // remove model immediately
}

double ModelLJShape::sigma(const int site_type, const ModelParams& params) {
  if (mixed_sigma_.size() == 0) {
    return params.select(sigma_index()).value(site_type);
  } else {
    return mixed_sigma_.value(site_type);
  }
}

double ModelLJShape::epsilon(const int site_type, const ModelParams& params) {
  if (mixed_epsilon_.size() == 0) {
    return params.select(epsilon_index()).value(site_type);
  } else {
    return mixed_epsilon_.value(site_type);
  }
}

double ModelLJShape::energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) {
  const double distance = -shape()->nearest_distance(wrapped_site);
  TRACE("distance: " << distance);
  const int type = site.type();
  const double cutoff = model_params.select(cutoff_index()).value(type);
  const double eps = epsilon(type, model_params);
  if (distance <= NEAR_ZERO && std::abs(eps) > NEAR_ZERO) {
    TRACE(MAX_PRECISION << distance << " " << eps);
    return NEAR_INFINITY;
  } else if (distance >= cutoff) {
    return 0.;
  } else {
    const double sig = sigma(type, model_params);
    const double en = energy(eps, sig, distance);
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

double ModelLJShapeEnergyAtCutoff::compute(const int type1, const ModelParams& model_params) {
  const double eps = model_->epsilon(type1, model_params);
  const double sig = model_->sigma(type1, model_params);
  const double cutoff = model_params.select("cutoff").value(type1);
  TRACE("eps " << MAX_PRECISION << eps);
  TRACE("sig " << sig);
  TRACE("cut " << cutoff);
  double en = 0.;
  if (cutoff > 0) {
    en = model_->energy(eps, sig, cutoff);
  }
  TRACE("en " << en);
  return en;
}

}  // namespace feasst
