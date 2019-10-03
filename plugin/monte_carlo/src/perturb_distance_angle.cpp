
#include "monte_carlo/include/perturb_distance_angle.h"

namespace feasst {

PerturbDistanceAngle::PerturbDistanceAngle(const argtype& args)
  : PerturbDistance(args) {
  class_name_ = "PerturbDistanceAngle";
}

class MapPerturbDistanceAngle {
 public:
  MapPerturbDistanceAngle() {
    auto obj = MakePerturbDistanceAngle();
    obj->deserialize_map()["PerturbDistanceAngle"] = obj;
  }
};

static MapPerturbDistanceAngle mapper_ = MapPerturbDistanceAngle();

std::shared_ptr<Perturb> PerturbDistanceAngle::create(std::istream& istr) const {
  return std::make_shared<PerturbDistanceAngle>(istr);
}

PerturbDistanceAngle::PerturbDistanceAngle(std::istream& istr)
  : PerturbDistance(istr) {
  ASSERT(class_name_ == "PerturbDistanceAngle", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(788 == version, "mismatch version: " << version);
  feasst_deserialize(&angle_, istr);
}

void PerturbDistanceAngle::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_distance_(ostr);
  feasst_serialize_version(788, ostr);
  feasst_serialize(angle_, ostr);
}

}  // namespace feasst
