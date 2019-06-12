#include "steppers/include/check.h"

namespace feasst {

class MapCheck {
 public:
  MapCheck() {
    Check().deserialize_map()["Check"] = MakeCheck();
  }
};

static MapCheck mapper_ = MapCheck();

}  // namespace feasst
