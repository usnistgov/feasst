#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/model_params.h"
#include "ewald/include/charge_screened.h"

namespace feasst {

class MapChargeScreened {
 public:
  MapChargeScreened() {
    ChargeScreened().deserialize_map()["ChargeScreened"] = MakeChargeScreened();
  }
};

static MapChargeScreened map_charge_screened_ = MapChargeScreened();

void ChargeScreened::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(4616, ostr);
  feasst_serialize(alpha_, ostr);
  feasst_serialize(conversion_factor_, ostr);
}

ChargeScreened::ChargeScreened(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 4616, "unrecognized verison: " << version);
  feasst_deserialize(&alpha_, istr);
  feasst_deserialize(&conversion_factor_, istr);
}

double ChargeScreened::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) const {
  const double distance = std::sqrt(squared_distance);
  if (std::abs(distance) < NEAR_ZERO) {
    return NEAR_INFINITY;
  }
  const double mixed_charge = model_params.mixed_charge()[type1][type2];
  return mixed_charge*conversion_factor_*erfc(alpha_*distance)/distance;
}

void ChargeScreened::precompute(const ModelParams& existing) {
  alpha_ = existing.property("alpha");
  conversion_factor_ = existing.constants().charge_conversion();
}

}  // namespace feasst
