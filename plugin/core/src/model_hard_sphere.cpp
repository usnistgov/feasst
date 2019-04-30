#include "core/include/model_hard_sphere.h"

namespace feasst {

class MapModelHardSphere {
 public:
  MapModelHardSphere() {
    ModelHardSphere().deserialize_map()["ModelHardSphere"] = std::make_shared<ModelHardSphere>();
  }
};

static MapModelHardSphere mapper_ = MapModelHardSphere();

}  // namespace feasst
