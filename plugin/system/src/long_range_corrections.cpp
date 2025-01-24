#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/select.h"
#include "configuration/include/model_params.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/long_range_corrections.h"

namespace feasst {

FEASST_MAPPER(LongRangeCorrections,);

void LongRangeCorrections::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
  feasst_serialize_version(875, ostr);
  feasst_serialize(mie_lambda_r_index_, ostr);
  feasst_serialize(mie_lambda_a_index_, ostr);
}

LongRangeCorrections::LongRangeCorrections(std::istream& istr)
  : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 874 && version <= 875, "version mismatch:" << version);
  if (version >= 875) {
    feasst_deserialize(&mie_lambda_r_index_, istr);
    feasst_deserialize(&mie_lambda_a_index_, istr);
  }
}

void LongRangeCorrections::precompute(Configuration * config) {
  VisitModel::precompute(config);
  ASSERT(config->domain().dimension() == 3, "LongRangeCorrections assumes 3D");
  mie_lambda_r_index_ = config->model_params().index("mie_lambda_r");
  mie_lambda_a_index_ = config->model_params().index("mie_lambda_a");
}

void LongRangeCorrections::compute(
    ModelOneBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  config->num_sites_of_type(group_index, &num_of_site_type_);
  config->num_sites_of_type(selection, &select_types_);
  DEBUG("sel: " << selection.str());
//  DEBUG("num sites of type in selection: " << feasst_str(select_types_));
  double en = 0.;
  double factor = -1;
  DEBUG("trial state " << selection.trial_state());
  if (selection.trial_state() == 3) {
    factor = 1.;
  }
  for (int type1 = 0; type1 < config->num_site_types(); ++type1) {
    DEBUG("type1 " << type1);
    const double num_type1 = num_of_site_type_[type1];
    const double num_type1_sel = select_types_[type1];
    DEBUG("num_particles " << config->num_particles());
    DEBUG("num_type1 " << num_type1);
    DEBUG("num_type1_sel " << num_type1_sel);
    for (int type2 = 0; type2 < config->num_site_types(); ++type2) {
      DEBUG("type2 " << type2);
      const double num_type2 = num_of_site_type_[type2];
      const double num_type2_sel = select_types_[type2];
      DEBUG("num_type2 " << num_type2);
      DEBUG("num_type2_sel " << num_type2_sel);
      en += (num_type1*num_type2_sel +
             num_type1_sel*num_type2 +
             factor*num_type1_sel*num_type2_sel)
        *energy_(type1, type2, config, model_params);
    }
  }
  set_energy(en);
}

void LongRangeCorrections::compute(
    ModelOneBody * model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  double en = 0;
  config->num_sites_of_type(group_index, &num_of_site_type_);
//  DEBUG("num sites of type in group: " << feasst_str(num_of_site_type_));
  for (int type1 = 0; type1 < config->num_site_types(); ++type1) {
    for (int type2 = 0; type2 < config->num_site_types(); ++type2) {
      en += num_of_site_type_[type1]*num_of_site_type_[type2]*
        energy_(type1, type2, config, model_params);
    }
  }
  TRACE("en " << en);
  set_energy(en);
}

double LongRangeCorrections::energy_(
    const int type1,
    const int type2,
    const Configuration * config,
    const ModelParams& model_params) const {
  DEBUG("type1: " << type1);
  DEBUG("type2: " << type2);
  const double epsilon = model_params.select(epsilon_index()).mixed_values()[type1][type2];
  DEBUG("epsilon: " << epsilon);
  const double sigma = model_params.select(sigma_index()).mixed_values()[type1][type2];
  DEBUG("sigma: " << sigma);
  const double cutoff = model_params.select(cutoff_index()).mixed_values()[type1][type2];
  DEBUG("cutoff: " << cutoff);
  if (std::abs(cutoff) < NEAR_ZERO || std::abs(epsilon) < NEAR_ZERO) {
    return 0;
  }
  if (mie_lambda_r_index_ == -1 || mie_lambda_a_index_ == -1) {
    const double prefactor = epsilon*(8./3.)*PI*std::pow(sigma, 3)*
      ((1./3.)*std::pow(sigma/cutoff, 9) - std::pow(sigma/cutoff, 3));
    DEBUG("prefactor: " << prefactor);
    const double en = prefactor/config->domain().volume();
    DEBUG("en: " << en);
    return en;
  }
  const double n = model_params.select(mie_lambda_r_index_).mixed_values()[type1][type2];
  const double m = model_params.select(mie_lambda_a_index_).mixed_values()[type1][type2];
  TRACE("n " << n << " m " << m);
  double prefactor = n/(n-m)*std::pow(n/m, m/(n-m));
  TRACE("prefactor: " << prefactor);
  prefactor *= 2.*epsilon*PI*std::pow(sigma, 3)/(m-3)*
    (((m-3.)/(n-3.))*std::pow(sigma/cutoff, n-3) - std::pow(sigma/cutoff, m-3));
  TRACE("prefactor: " << prefactor);
  const double en = prefactor/config->domain().volume();
  TRACE("en: " << en);
  return en;
}

}  // namespace feasst
