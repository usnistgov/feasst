#include "models/include/square_well.h"

namespace feasst {

class MapSquareWell {
 public:
  MapSquareWell() {
    SquareWell().deserialize_map()["SquareWell"] = MakeSquareWell();
  }
};

static MapSquareWell mapper_ = MapSquareWell();

}  // namespace feasst
