#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/model_params.h"
#include "charge/include/debye_huckel.h"

namespace feasst {

class MapDebyeHuckel {
 public:
  MapDebyeHuckel() {
    auto obj = MakeDebyeHuckel({{"kappa", "1"}, {"dielectric", "1"}});
    obj->deserialize_map()["DebyeHuckel"] = obj;
  }
};

static MapDebyeHuckel map_charge_screened_ = MapDebyeHuckel();

DebyeHuckel::DebyeHuckel(argtype * args) {
  class_name_ = "DebyeHuckel";
  kappa_ = dble("kappa", args);
  dielectric_ = dble("dielectric", args);
}
DebyeHuckel::DebyeHuckel(argtype args) : DebyeHuckel(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void DebyeHuckel::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(3682, ostr);
  feasst_serialize(conversion_factor_, ostr);
  feasst_serialize(kappa_, ostr);
  feasst_serialize(dielectric_, ostr);
}

DebyeHuckel::DebyeHuckel(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3682, "unrecognized verison: " << version);
  feasst_deserialize(&conversion_factor_, istr);
  feasst_deserialize(&kappa_, istr);
  feasst_deserialize(&dielectric_, istr);
}

double DebyeHuckel::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  const double distance = std::sqrt(squared_distance);
  TRACE("distance " << distance);
  if (std::abs(distance) < NEAR_ZERO) {
    TRACE("near inf");
    return NEAR_INFINITY;
  }
  const double mixed_charge = model_params.select(charge_index()).mixed_values()[type1][type2];
  //TRACE("mixed_charge " << mixed_charge);
  //TRACE("conversion_factor_ " << conversion_factor_);
  return mixed_charge*conversion_factor_*std::exp(-kappa_*distance)/distance/dielectric_;
}

void DebyeHuckel::precompute(const ModelParams& existing) {
  Model::precompute(existing);
  conversion_factor_ = existing.constants().charge_conversion();
}

}  // namespace feasst
