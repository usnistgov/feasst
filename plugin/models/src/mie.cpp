#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/configuration.h"
#include "models/include/mie.h"

namespace feasst {

Mie::Mie(argtype * args) {
  class_name_ = "Mie";
}
Mie::Mie(argtype args) : Mie(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(Mie,);

void Mie::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_mie_(ostr);
}

void Mie::serialize_mie_(std::ostream& ostr) const {
  serialize_model_(ostr);
  feasst_serialize_version(2905, ostr);
  feasst_serialize(mie_lambda_r_index_, ostr);
  feasst_serialize(mie_lambda_a_index_, ostr);
  feasst_serialize_fstobj(prefactor_, ostr);
}

Mie::Mie(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(2905 == version, version);
  feasst_deserialize(&mie_lambda_r_index_, istr);
  feasst_deserialize(&mie_lambda_a_index_, istr);
  feasst_deserialize_fstobj(&prefactor_, istr);
}

void Mie::precompute(const Configuration& config) {
  Model::precompute(config);
  const ModelParams& existing = config.model_params();
  mie_lambda_r_index_ = existing.index("mie_lambda_r");
  mie_lambda_a_index_ = existing.index("mie_lambda_a");
  ASSERT(mie_lambda_r_index_ != -1 && mie_lambda_a_index_ != -1,
    "Mie potential requires Site Properties mie_lambda_r and mie_lambda_a");
  prefactor_.set_param(existing);
  for (int type1 = 0; type1 < existing.size(); ++type1) {
    for (int type2 = 0; type2 < existing.size(); ++type2) {
      prefactor_.compute(type1, type2, existing);
    }
  }
}

double Mie::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  const double sigma = model_params.select(sigma_index()).mixed_values()[type1][type2];
  const double s_r_sq = sigma*sigma/squared_distance;
  const double n = model_params.select(mie_lambda_r_index_).mixed_values()[type1][type2];
  TRACE("n " << n);
  const double m = model_params.select(mie_lambda_a_index_).mixed_values()[type1][type2];
  TRACE("m " << m);
  const double prefactor = prefactor_.mixed_value(type1, type2);
  TRACE("prefactor " << prefactor);
  const double en = prefactor*(std::pow(s_r_sq, 0.5*n) - std::pow(s_r_sq, 0.5*m));
  TRACE("en " << en);
  return en;
}

}  // namespace feasst
