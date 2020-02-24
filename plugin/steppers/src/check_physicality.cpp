#include "steppers/include/check_physicality.h"

namespace feasst {

class MapCheckPhysicality {
 public:
  MapCheckPhysicality() {
    CheckPhysicality().deserialize_map()["CheckPhysicality"] =
      MakeCheckPhysicality();
  }
};

static MapCheckPhysicality mapper_ = MapCheckPhysicality();

}  // namespace feasst
