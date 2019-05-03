#include "models/include/model_square_well.h"

namespace feasst {

class MapModelSquareWell {
 public:
  MapModelSquareWell() {
    ModelSquareWell().deserialize_map()["ModelSquareWell"] = MakeModelSquareWell();
  }
};

static MapModelSquareWell mapper_ = MapModelSquareWell();

}  // namespace feasst
