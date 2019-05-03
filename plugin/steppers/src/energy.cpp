#include "steppers/include/energy.h"

namespace feasst {

Energy::Energy(const argtype &args) : Analyze(args) {}

class MapEnergy {
 public:
  MapEnergy() {
    auto obj = MakeEnergy();
    obj->deserialize_map()["Energy"] = obj;
  }
};

static MapEnergy mapper_ = MapEnergy();

}  // namespace feasst
