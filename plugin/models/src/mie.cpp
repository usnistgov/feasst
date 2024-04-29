#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "models/include/mie.h"

namespace feasst {

Mie::Mie(argtype * args) {
  class_name_ = "Mie";
}
Mie::Mie(argtype args) : Mie(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapMie {
 public:
  MapMie() {
    Mie().deserialize_map()["Mie"] = MakeMie();
  }
};

static MapMie map_mie_ = MapMie();

void Mie::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_mie_(ostr);
}

void Mie::serialize_mie_(std::ostream& ostr) const {
  serialize_model_(ostr);
  feasst_serialize_version(2905, ostr);
  feasst_serialize(mie_lambda_r_index_, ostr);
  feasst_serialize(mie_lambda_a_index_, ostr);
}

Mie::Mie(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(2905 == version, version);
  feasst_deserialize(&mie_lambda_r_index_, istr);
  feasst_deserialize(&mie_lambda_a_index_, istr);
}

void Mie::precompute(const ModelParams& existing) {
  Model::precompute(existing);
  mie_lambda_r_index_ = existing.index("mie_lambda_r");
  mie_lambda_a_index_ = existing.index("mie_lambda_a");
  ASSERT(mie_lambda_r_index_ != -1 && mie_lambda_a_index_ != -1,
    "Mie potential requires Site Properties mie_lambda_r and mie_lambda_a");
}

double Mie::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  const double sigma = model_params.select(sigma_index()).mixed_values()[type1][type2];
  const double epsilon = model_params.select(epsilon_index()).mixed_values()[type1][type2];
  const double s_r = sigma/std::sqrt(squared_distance);
  const double n = model_params.select(mie_lambda_r_index_).mixed_values()[type1][type2];
  const double m = model_params.select(mie_lambda_a_index_).mixed_values()[type1][type2];
  TRACE("n " << n);
  TRACE("m " << m);
  const double prefactor = n/(n-m)*std::pow(n/m, m/(n-m));
  const double en = prefactor*epsilon*(std::pow(s_r, n) - std::pow(s_r, m));
  return en;
}

}  // namespace feasst
