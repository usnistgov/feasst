
#include "monte_carlo/include/perturb_add.h"

namespace feasst {

PerturbAdd::PerturbAdd(const argtype& args) : Perturb(args) {
  class_name_ = "PerturbAdd";
  disable_tunable_();
}

class MapPerturbAdd {
 public:
  MapPerturbAdd() {
    auto obj = MakePerturbAdd();
    obj->deserialize_map()["PerturbAdd"] = obj;
  }
};

static MapPerturbAdd mapper_ = MapPerturbAdd();

std::shared_ptr<Perturb> PerturbAdd::create(std::istream& istr) const {
  return std::make_shared<PerturbAdd>(istr);
}

PerturbAdd::PerturbAdd(std::istream& istr)
  : Perturb(istr) {
  ASSERT(class_name_ == "PerturbAdd", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(730 == version, "mismatch version: " << version);
}

void PerturbAdd::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(730, ostr);
}

}  // namespace feasst
