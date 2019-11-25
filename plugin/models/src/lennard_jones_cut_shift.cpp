#include "models/include/lennard_jones_cut_shift.h"

namespace feasst {

class MapLennardJonesCutShift {
 public:
  MapLennardJonesCutShift() {
    LennardJonesCutShift().deserialize_map()["LennardJonesCutShift"] = MakeLennardJonesCutShift();
  }
};

static MapLennardJonesCutShift map_lennard_jones_cut_shift_ = MapLennardJonesCutShift();

}  // namespace feasst
