#include "core/include/model_lj.h"

namespace feasst {

class MapModelLJ {
 public:
  MapModelLJ() {
    ModelLJ().deserialize_map()["ModelLJ"] = MakeModelLJ();
  }
};

static MapModelLJ mapper_ = MapModelLJ();

void ModelLJ::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_lj_(ostr);
}

void ModelLJ::serialize_model_lj_(std::ostream& ostr) const {
  feasst_serialize_version(763, ostr);
  feasst_serialize(hard_sphere_threshold_sq_, ostr);
}

ModelLJ::ModelLJ(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(763 == version, version);
  feasst_deserialize(&hard_sphere_threshold_sq_, istr);
}

}  // namespace feasst
