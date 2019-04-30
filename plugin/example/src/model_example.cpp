#include "example/include/model_example.h"

namespace feasst {

class MapModelExample {
 public:
  MapModelExample() {
    ModelExample().deserialize_map()["ModelExample"] = MakeModelExample();
  }
};

static MapModelExample mapper_ = MapModelExample();

ModelExample::ModelExample() {}

}  // namespace feasst
