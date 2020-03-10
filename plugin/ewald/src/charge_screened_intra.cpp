#include "ewald/include/charge_screened_intra.h"

namespace feasst {

class MapChargeScreenedIntra {
 public:
  MapChargeScreenedIntra() {
    ChargeScreenedIntra().deserialize_map()["ChargeScreenedIntra"] = MakeChargeScreenedIntra();
  }
};

static MapChargeScreenedIntra map_charge_screened_intra_ = MapChargeScreenedIntra();

void ChargeScreenedIntra::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(304, ostr);
  feasst_serialize(alpha_, ostr);
  feasst_serialize(conversion_factor_, ostr);
}

ChargeScreenedIntra::ChargeScreenedIntra(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 304, "unrecognized verison: " << version);
  feasst_deserialize(&alpha_, istr);
  feasst_deserialize(&conversion_factor_, istr);
}

}  // namespace feasst
