#include "ewald/include/model_charge_intra.h"

namespace feasst {

class MapModelChargeIntra {
 public:
  MapModelChargeIntra() {
    ModelChargeIntra().deserialize_map()["ModelChargeIntra"] = MakeModelChargeIntra();
  }
};

static MapModelChargeIntra map_model_charge_intra_ = MapModelChargeIntra();

}  // namespace feasst
