#include "steppers/include/energy.h"

namespace feasst {

class MapEnergy {
 public:
  MapEnergy() {
    auto obj = MakeEnergy();
    obj->deserialize_map()["Energy"] = obj;
  }
};

static MapEnergy mapper_ = MapEnergy();

Energy::Energy(const argtype &args) : Analyze(args) {
  args_.init(args);
  energy_.set_block(args_.key("num_blocks").dflt(str(1e5)).integer());
}

}  // namespace feasst
