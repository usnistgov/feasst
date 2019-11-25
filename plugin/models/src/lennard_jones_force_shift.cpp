#include "models/include/lennard_jones_force_shift.h"

namespace feasst {

class MapLennardJonesForceShift {
 public:
  MapLennardJonesForceShift() {
    LennardJonesForceShift().deserialize_map()["LennardJonesForceShift"] = MakeLennardJonesForceShift();
  }
};

static MapLennardJonesForceShift map_lennard_jones_force_shift_ = MapLennardJonesForceShift();

}  // namespace feasst
