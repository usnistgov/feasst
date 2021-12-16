#include <cmath>
#include "utils/include/serialize.h"
#include "configuration/include/model_params.h"
#include "models/include/yukawa.h"

namespace feasst {

class MapYukawa {
 public:
  MapYukawa() {
    Yukawa().deserialize_map()["Yukawa"] = MakeYukawa();
  }
};

static MapYukawa map_model_lj_alpha_ = MapYukawa();

Yukawa::Yukawa() {
  class_name_ = "Yukawa";
  set_kappa();
}

void Yukawa::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(6505, ostr);
  feasst_serialize(kappa_, ostr);
}

Yukawa::Yukawa(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6505, "unrecognized version: " << version);
  feasst_deserialize(&kappa_, istr);
}

double Yukawa::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  const double epsilon = model_params.select(epsilon_index()).mixed_values()[type1][type2];
  const double sigma = model_params.select(sigma_index()).mixed_values()[type1][type2];
  const double distance = sqrt(squared_distance);
  return epsilon*std::exp(-kappa_*(distance/sigma - 1.))/(distance/sigma);
}

}  // namespace feasst
