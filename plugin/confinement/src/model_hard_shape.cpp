#include "confinement/include/model_hard_shape.h"

namespace feasst {

class MapModelHardShape {
 public:
  MapModelHardShape() {
    ModelHardShape().deserialize_map()["ModelHardShape"] = MakeModelHardShape();
  }
};

static MapModelHardShape map_model_hard_shape_ = MapModelHardShape();

}  // namespace feasst
