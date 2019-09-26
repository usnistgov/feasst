
#include "monte_carlo/include/perturb_remove.h"

namespace feasst {

class MapPerturbRemove {
 public:
  MapPerturbRemove() {
    auto obj = MakePerturbRemove();
    obj->deserialize_map()["PerturbRemove"] = obj;
  }
};

static MapPerturbRemove mapper_ = MapPerturbRemove();

std::shared_ptr<Perturb> PerturbRemove::create(std::istream& istr) const {
  return std::make_shared<PerturbRemove>(istr);
}

PerturbRemove::PerturbRemove(std::istream& istr)
  : Perturb(istr) {
  ASSERT(class_name_ == "PerturbRemove", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(143 == version, "mismatch version: " << version);
}

void PerturbRemove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(143, ostr);
}

}  // namespace feasst
