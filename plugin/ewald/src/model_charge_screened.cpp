#include "ewald/include/model_charge_screened.h"

namespace feasst {

class MapModelChargeScreened {
 public:
  MapModelChargeScreened() {
    ModelChargeScreened().deserialize_map()["ModelChargeScreened"] = MakeModelChargeScreened();
  }
};

static MapModelChargeScreened map_model_charge_screened_ = MapModelChargeScreened();

}  // namespace feasst
