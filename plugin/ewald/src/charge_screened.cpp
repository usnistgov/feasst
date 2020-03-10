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

}  // namespace feasst
