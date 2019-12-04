#include "ewald/include/charge_screened_intra.h"

namespace feasst {

class MapChargeScreenedIntra {
 public:
  MapChargeScreenedIntra() {
    ChargeScreenedIntra().deserialize_map()["ChargeScreenedIntra"] = MakeChargeScreenedIntra();
  }
};

static MapChargeScreenedIntra map_charge_screened_intra_ = MapChargeScreenedIntra();

}  // namespace feasst
