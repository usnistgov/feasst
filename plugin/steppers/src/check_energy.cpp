#include "steppers/include/check_energy.h"

namespace feasst {

class MapCheckEnergy {
 public:
  MapCheckEnergy() {
    CheckEnergy().deserialize_map()["CheckEnergy"] = MakeCheckEnergy();
  }
};

static MapCheckEnergy mapper_energy_check_ = MapCheckEnergy();

}  // namespace feasst
