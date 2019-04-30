#include "models/include/model_lj_alpha.h"

namespace feasst {

ModelLJAlpha::ModelLJAlpha() {
  set_alpha();
}

class MapModelLJAlpha {
 public:
  MapModelLJAlpha() {
    ModelLJAlpha().deserialize_map()["ModelLJAlpha"] = MakeModelLJAlpha();
  }
};

static MapModelLJAlpha map_model_lj_alpha_ = MapModelLJAlpha();

}  // namespace feasst
