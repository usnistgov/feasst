#include "system/include/model_empty.h"

namespace feasst {

class MapModelEmpty {
 public:
  MapModelEmpty() {
    ModelEmpty().deserialize_map()["ModelEmpty"]
      = std::make_shared<ModelEmpty>();
  }
};

static MapModelEmpty mapper_ = MapModelEmpty();

}  // namespace feasst
