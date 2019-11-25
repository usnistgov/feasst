#include "models/include/lennard_jones_alpha.h"

namespace feasst {

LennardJonesAlpha::LennardJonesAlpha() {
  set_alpha();
}

class MapLennardJonesAlpha {
 public:
  MapLennardJonesAlpha() {
    LennardJonesAlpha().deserialize_map()["LennardJonesAlpha"] = MakeLennardJonesAlpha();
  }
};

static MapLennardJonesAlpha map_lennard_jones_alpha_ = MapLennardJonesAlpha();

}  // namespace feasst
