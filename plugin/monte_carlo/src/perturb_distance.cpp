
#include "monte_carlo/include/perturb_distance.h"

namespace feasst {

class MapPerturbDistance {
 public:
  MapPerturbDistance() {
    auto obj = MakePerturbDistance();
    obj->deserialize_map()["PerturbDistance"] = obj;
  }
};

static MapPerturbDistance mapper_ = MapPerturbDistance();

std::shared_ptr<Perturb> PerturbDistance::create(std::istream& istr) const {
  return std::make_shared<PerturbDistance>(istr);
}

PerturbDistance::PerturbDistance(std::istream& istr)
  : PerturbMove(istr) {
  // HWH can't check this if this is a base class
  // ASSERT(class_name_ == "PerturbDistance", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(228 == version, "mismatch version: " << version);
  feasst_deserialize(&distance_, istr);
}

void PerturbDistance::serialize_perturb_distance_(std::ostream& ostr) const {
  serialize_perturb_(ostr);
  feasst_serialize_version(228, ostr);
  feasst_serialize(distance_, ostr);
}

void PerturbDistance::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_distance_(ostr);
}

}  // namespace feasst
