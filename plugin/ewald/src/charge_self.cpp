#include "ewald/include/charge_self.h"

namespace feasst {

class MapChargeSelf {
 public:
  MapChargeSelf() {
    ChargeSelf().deserialize_map()["ChargeSelf"] = MakeChargeSelf();
  }
};

static MapChargeSelf map_charge_self_ = MapChargeSelf();

void ChargeSelf::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(8833, ostr);
  feasst_serialize(alpha_, ostr);
  feasst_serialize(conversion_factor_, ostr);
}

ChargeSelf::ChargeSelf(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8833, "unrecognized verison: " << version);
  feasst_deserialize(&alpha_, istr);
  feasst_deserialize(&conversion_factor_, istr);
}

}  // namespace feasst
