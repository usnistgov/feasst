#include "confinement/include/model_hard_shape.h"

namespace feasst {

class MapModelHardShape {
 public:
  MapModelHardShape() {
    ModelHardShape().deserialize_map()["ModelHardShape"] =
      std::make_shared<ModelHardShape>();
  }
};

static MapModelHardShape map_model_hard_shape_ = MapModelHardShape();

}  // namespace feasst
