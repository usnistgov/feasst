#include "models/include/model_lj_cut_shift.h"

namespace feasst {

class MapModelLJCutShift {
 public:
  MapModelLJCutShift() {
    ModelLJCutShift().deserialize_map()["ModelLJCutShift"] = MakeModelLJCutShift();
  }
};

static MapModelLJCutShift map_model_lj_cut_shift_ = MapModelLJCutShift();

}  // namespace feasst
