
#include "monte_carlo/include/perturb_anywhere.h"

namespace feasst {

class MapPerturbAnywhere {
 public:
  MapPerturbAnywhere() {
    auto obj = MakePerturbAnywhere();
    obj->deserialize_map()["PerturbAnywhere"] = obj;
  }
};

static MapPerturbAnywhere mapper_ = MapPerturbAnywhere();

std::shared_ptr<Perturb> PerturbAnywhere::create(std::istream& istr) const {
  return std::make_shared<PerturbAnywhere>(istr);
}

PerturbAnywhere::PerturbAnywhere(std::istream& istr)
  : PerturbMove(istr) {
  ASSERT(class_name_ == "PerturbAnywhere", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(461 == version, "mismatch version: " << version);
  feasst_deserialize_fstobj(&rotate_, istr);
}

void PerturbAnywhere::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(461, ostr);
  feasst_serialize_fstobj(rotate_, ostr);
}

}  // namespace feasst
