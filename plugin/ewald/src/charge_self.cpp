#include "ewald/include/charge_self.h"

namespace feasst {

class MapChargeSelf {
 public:
  MapChargeSelf() {
    ChargeSelf().deserialize_map()["ChargeSelf"] = MakeChargeSelf();
  }
};

static MapChargeSelf map_charge_self_ = MapChargeSelf();

}  // namespace feasst
