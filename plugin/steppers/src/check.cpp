#include "steppers/include/check.h"

namespace feasst {

class MapCheck {
 public:
  MapCheck() {
    Check().deserialize_map()["Check"] = MakeCheck();
  }
};

static MapCheck mapper_ = MapCheck();

class MapEnergyCheck {
 public:
  MapEnergyCheck() {
    EnergyCheck().deserialize_map()["EnergyCheck"] = MakeEnergyCheck();
  }
};

static MapEnergyCheck mapper_energy_check_ = MapEnergyCheck();

}  // namespace feasst
