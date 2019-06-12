#include "system/include/model_two_body_factory.h"

namespace feasst {

class MapModelTwoBodyFactory {
 public:
  MapModelTwoBodyFactory() {
    ModelTwoBodyFactory().deserialize_map()["ModelTwoBodyFactory"] =
      MakeModelTwoBodyFactory();
  }
};

static MapModelTwoBodyFactory mapper_ = MapModelTwoBodyFactory();

}  // namespace feasst
