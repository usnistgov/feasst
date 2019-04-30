#include "ewald/include/model_charge_self.h"

namespace feasst {

class MapModelChargeSelf {
 public:
  MapModelChargeSelf() {
    ModelChargeSelf().deserialize_map()["ModelChargeSelf"] = MakeModelChargeSelf();
  }
};

static MapModelChargeSelf map_model_charge_self_ = MapModelChargeSelf();

}  // namespace feasst
