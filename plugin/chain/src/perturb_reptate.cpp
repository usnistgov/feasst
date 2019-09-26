
#include "chain/include/perturb_reptate.h"

namespace feasst {

class MapPerturbReptate {
 public:
  MapPerturbReptate() {
    auto obj = MakePerturbReptate();
    obj->deserialize_map()["PerturbReptate"] = obj;
  }
};

static MapPerturbReptate mapper_ = MapPerturbReptate();

std::shared_ptr<Perturb> PerturbReptate::create(std::istream& istr) const {
  return std::make_shared<PerturbReptate>(istr);
}

PerturbReptate::PerturbReptate(std::istream& istr)
  : PerturbDistance(istr) {
  ASSERT(class_name_ == "PerturbReptate", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(863 == version, "mismatch version: " << version);
}

void PerturbReptate::serialize_perturb_reptate_(std::ostream& ostr) const {
  serialize_perturb_distance_(ostr);
  feasst_serialize_version(863, ostr);
}

void PerturbReptate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_reptate_(ostr);
}

}  // namespace feasst
