#include "models/include/model_lj_force_shift.h"

namespace feasst {

class MapModelLJForceShift {
 public:
  MapModelLJForceShift() {
    ModelLJForceShift().deserialize_map()["ModelLJForceShift"] = MakeModelLJForceShift();
  }
};

static MapModelLJForceShift map_model_lj_force_shift_ = MapModelLJForceShift();

}  // namespace feasst
