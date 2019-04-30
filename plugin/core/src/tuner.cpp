#include "core/include/tuner.h"

namespace feasst {

class MapTuner {
 public:
  MapTuner() {
    Tuner().deserialize_map()["Tuner"] = MakeTuner();
  }
};

static MapTuner mapper_ = MapTuner();

}  // namespace feasst
