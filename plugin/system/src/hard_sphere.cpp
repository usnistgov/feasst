#include "system/include/hard_sphere.h"

namespace feasst {

class MapHardSphere {
 public:
  MapHardSphere() {
    HardSphere().deserialize_map()["HardSphere"] = std::make_shared<HardSphere>();
  }
};

static MapHardSphere mapper_ = MapHardSphere();

}  // namespace feasst
