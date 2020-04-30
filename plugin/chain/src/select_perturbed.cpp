#include "chain/include/select_perturbed.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapSelectPerturbed {
 public:
  MapSelectPerturbed() {
    auto obj = MakeSelectPerturbed();
    obj->deserialize_map()["SelectPerturbed"] = obj;
  }
};

static MapSelectPerturbed mapper_ = MapSelectPerturbed();

std::shared_ptr<TrialSelect> SelectPerturbed::create(std::istream& istr) const {
  return std::make_shared<SelectPerturbed>(istr);
}

SelectPerturbed::SelectPerturbed(std::istream& istr)
  : TrialSelect(istr) {
  // ASSERT(class_name_ == "SelectPerturbed", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(607 == version, "mismatch version: " << version);
}

void SelectPerturbed::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_(ostr);
  feasst_serialize_version(607, ostr);
}

}  // namespace feasst
