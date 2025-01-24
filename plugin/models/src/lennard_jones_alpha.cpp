#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "models/include/lennard_jones_alpha.h"

namespace feasst {

LennardJonesAlpha::LennardJonesAlpha(argtype args) : LennardJonesAlpha(&args) {
  feasst_check_all_used(args);
}
LennardJonesAlpha::LennardJonesAlpha(argtype * args) : LennardJones(args) {
  class_name_ = "LennardJonesAlpha";
  alpha_ = dble("alpha", args, 6);
  two_raised_inv_alpha_ = std::pow(2, 1./alpha_);
  DEBUG("two_raised_inv_alpha_ " << two_raised_inv_alpha_);
}

FEASST_MAPPER(LennardJonesAlpha,);

void LennardJonesAlpha::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_lennard_jones_alpha_(ostr);
}

void LennardJonesAlpha::serialize_lennard_jones_alpha_(
    std::ostream& ostr) const {
  serialize_lennard_jones_(ostr);
  feasst_serialize_version(714, ostr);
  feasst_serialize(alpha_, ostr);
  feasst_serialize(delta_sigma_index_, ostr);
  feasst_serialize(lambda_index_, ostr);
  feasst_serialize(two_raised_inv_alpha_, ostr);
}

LennardJonesAlpha::LennardJonesAlpha(std::istream& istr) : LennardJones(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 713 && version <= 714, version);
  feasst_deserialize(&alpha_, istr);
  feasst_deserialize(&delta_sigma_index_, istr);
  if (version < 714) {
    double lambda;
    feasst_deserialize(&lambda, istr);
  }
  feasst_deserialize(&lambda_index_, istr);
  feasst_deserialize(&two_raised_inv_alpha_, istr);
}

void LennardJonesAlpha::precompute(const ModelParams& existing) {
  Model::precompute(existing);
  delta_sigma_index_ = existing.index("delta_sigma");
  lambda_index_ = existing.index("lambda");
  DEBUG("lambda_index_ " << lambda_index_);
}

double LennardJonesAlpha::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  const double sigma = model_params.select(sigma_index()).mixed_values()[type1][type2];
  const double sigma_squared = sigma*sigma;
  TRACE("squared_distance " << squared_distance);
  TRACE("hard_sphere_threshold_sq " << hard_sphere_threshold_sq());
  TRACE("sigma_squared " << sigma_squared);
  if (squared_distance == 0 ||
      squared_distance < hard_sphere_threshold_sq()*sigma_squared) {
    return NEAR_INFINITY;
  }
  const double epsilon = model_params.select(epsilon_index()).mixed_values()[type1][type2];
  const double distance = std::sqrt(squared_distance);
  TRACE("distance " << distance);
  double rinv2;
  double delta_sigma = 0;
  if (delta_sigma_index_ == -1) {
    rinv2 = sigma_squared/squared_distance;
  } else {
    delta_sigma = model_params.select(delta_sigma_index_).mixed_values()[type1][type2];
    const double rinv = (sigma + delta_sigma)/(distance + delta_sigma);
    rinv2 = rinv*rinv;
  }
  const double rinv_alpha = std::pow(rinv2, 0.5*alpha_);
  const double en = 4.*epsilon*rinv_alpha*(rinv_alpha - 1.);
  TRACE("lambda_index? " << lambda_index_);
  if (lambda_index_ == -1) {
    return en;
  } else {
    const double lambda = model_params.select(lambda_index_).mixed_values()[type1][type2];
    TRACE("lambda " << lambda);
    const double rmin = two_raised_inv_alpha_*(delta_sigma + sigma) - delta_sigma;
    if (distance <= rmin) {
      return en + epsilon*(1 - lambda);
    } else {
      return lambda*en;
    }
  }
}

double LennardJonesAlpha::du_dr(
    const double distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) const {
  const double epsilon = model_params.select(epsilon_index()).mixed_values()[type1][type2];
  const double sigma = model_params.select(sigma_index()).mixed_values()[type1][type2];
  double delta_sigma = 0.;
  if (delta_sigma_index_ != -1) {
    delta_sigma = model_params.select(delta_sigma_index_).mixed_values()[type1][type2];
  }
  const double rinv = (sigma + delta_sigma)/(distance + delta_sigma);
  if (sigma == 0) {
    return 0.;
  }
  const double dudr = 4.*epsilon*alpha()/(sigma + delta_sigma)*
    (-2*std::pow(rinv, 2*alpha() + 1)
      + std::pow(rinv, alpha() + 1));
  const double rmin = two_raised_inv_alpha_*(delta_sigma + sigma) - delta_sigma;
  if (lambda_index_ == -1 || distance <= rmin) {
    return dudr;
  } else {
    const double lambda = model_params.select(lambda_index_).mixed_values()[type1][type2];
    return lambda*dudr;
  }
}

void LennardJonesAlpha::set_wca(const int site_type1, const int site_type2,
    ModelParams * params) const {
  const double sigma = params->select("sigma").mixed_values()[site_type1][site_type2];
  const double r_wca = std::pow(2, 1./alpha())*sigma;
  params->set("cutoff", site_type1, site_type2, r_wca);
}

double EnergyAtCutOff::compute(const int type1, const int type2,
    const ModelParams& model_params) {
  const double cutoff = model_params.select("cutoff").mixed_values()[type1][type2];
  return model_->energy_without_shift(cutoff*cutoff, type1, type2, model_params);
}

double EnergyDerivAtCutOff::compute(const int type1, const int type2,
    const ModelParams& model_params) {
  const double cutoff = model_params.select("cutoff").mixed_values()[type1][type2];
  return model_->du_dr(cutoff, type1, type2, model_params);
}

double Lambda::mix_(const double value1, const double value2) {
  return std::sqrt(value1*value2);
}

FEASST_MAPPER_RENAME(DeltaSigma, delta_sigma,);

void DeltaSigma::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(8573, ostr);
}

DeltaSigma::DeltaSigma(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8573, "mismatch version: " << version);
}

FEASST_MAPPER_RENAME(Lambda, lambda,);

void Lambda::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(2045, ostr);
}

Lambda::Lambda(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2045, "mismatch version: " << version);
}

}  // namespace feasst
