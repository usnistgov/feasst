
#include "monte_carlo/include/perturb_anywhere.h"

namespace feasst {

class MapPerturbTranslate {
 public:
  MapPerturbTranslate() {
    auto obj = MakePerturbTranslate();
    obj->deserialize_map()["PerturbTranslate"] = obj;
  }
};

static MapPerturbTranslate mapper_ = MapPerturbTranslate();

std::shared_ptr<Perturb> PerturbTranslate::create(std::istream& istr) const {
  return std::make_shared<PerturbTranslate>(istr);
}

PerturbTranslate::PerturbTranslate(std::istream& istr)
  : PerturbMove(istr) {
  ASSERT(class_name_ == "PerturbTranslate", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(564 == version, "mismatch version: " << version);
}

void PerturbTranslate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(564, ostr);
}

}  // namespace feasst
