#include "models/include/model_yukawa.h"

namespace feasst {

class MapModelYukawa {
 public:
  MapModelYukawa() {
    ModelYukawa().deserialize_map()["ModelYukawa"] = MakeModelYukawa();
  }
};

static MapModelYukawa map_model_lj_alpha_ = MapModelYukawa();

ModelYukawa::ModelYukawa() {
  set_kappa();
}

}  // namespace feasst
