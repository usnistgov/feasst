#include "core/include/wall_clock_limit.h"

namespace feasst {

// this example shows how to handle a required argument.
class MapWallClockLimit {
 public:
  MapWallClockLimit() {
    auto obj = MakeWallClockLimit({{"max_hours", "0"}});
    obj->deserialize_map()["WallClockLimit"] = obj;
  }
};

static MapWallClockLimit mapper_ = MapWallClockLimit();

}  // namespace feasst
