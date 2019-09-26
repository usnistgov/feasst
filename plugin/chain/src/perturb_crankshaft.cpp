
#include "chain/include/perturb_crankshaft.h"

namespace feasst {

class MapPerturbCrankshaft {
 public:
  MapPerturbCrankshaft() {
    auto obj = MakePerturbCrankshaft();
    obj->deserialize_map()["PerturbCrankshaft"] = obj;
  }
};

static MapPerturbCrankshaft mapper_ = MapPerturbCrankshaft();

std::shared_ptr<Perturb> PerturbCrankshaft::create(std::istream& istr) const {
  return std::make_shared<PerturbCrankshaft>(istr);
}

PerturbCrankshaft::PerturbCrankshaft(std::istream& istr)
  : PerturbRotate(istr) {
  ASSERT(class_name_ == "PerturbCrankshaft", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(555 == version, "mismatch version: " << version);
}

void PerturbCrankshaft::serialize_perturb_crankshaft_(std::ostream& ostr) const {
  serialize_perturb_rotate_(ostr);
  feasst_serialize_version(555, ostr);
}

void PerturbCrankshaft::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_crankshaft_(ostr);
}

}  // namespace feasst
