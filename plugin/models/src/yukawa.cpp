#include "models/include/yukawa.h"

namespace feasst {

class MapYukawa {
 public:
  MapYukawa() {
    Yukawa().deserialize_map()["Yukawa"] = MakeYukawa();
  }
};

static MapYukawa map_model_lj_alpha_ = MapYukawa();

Yukawa::Yukawa() {
  set_kappa();
}

}  // namespace feasst
