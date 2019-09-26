
#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

class MapPerturbRotate {
 public:
  MapPerturbRotate() {
    auto obj = MakePerturbRotate();
    obj->deserialize_map()["PerturbRotate"] = obj;
  }
};

static MapPerturbRotate mapper_ = MapPerturbRotate();

std::shared_ptr<Perturb> PerturbRotate::create(std::istream& istr) const {
  return std::make_shared<PerturbRotate>(istr);
}

PerturbRotate::PerturbRotate(std::istream& istr)
  : PerturbMove(istr) {
  // ASSERT(class_name_ == "PerturbRotate", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(448 == version, "mismatch version: " << version);
}

void PerturbRotate::serialize_perturb_rotate_(std::ostream& ostr) const {
  serialize_perturb_(ostr);
  feasst_serialize_version(448, ostr);
}

void PerturbRotate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_rotate_(ostr);
}

}  // namespace feasst
