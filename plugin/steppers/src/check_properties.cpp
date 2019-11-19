#include "steppers/include/check_properties.h"

namespace feasst {

class MapCheckProperties {
 public:
  MapCheckProperties() {
    CheckProperties().deserialize_map()["CheckProperties"] = MakeCheckProperties();
  }
};

static MapCheckProperties mapper_check_properties_ = MapCheckProperties();

}  // namespace feasst
