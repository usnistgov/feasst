#include "chain/include/trial_select_perturbed.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapTrialSelectPerturbed {
 public:
  MapTrialSelectPerturbed() {
    auto obj = MakeTrialSelectPerturbed();
    obj->deserialize_map()["TrialSelectPerturbed"] = obj;
  }
};

static MapTrialSelectPerturbed mapper_ = MapTrialSelectPerturbed();

std::shared_ptr<TrialSelect> TrialSelectPerturbed::create(std::istream& istr) const {
  return std::make_shared<TrialSelectPerturbed>(istr);
}

TrialSelectPerturbed::TrialSelectPerturbed(std::istream& istr)
  : TrialSelect(istr) {
  // ASSERT(class_name_ == "TrialSelectPerturbed", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(607 == version, "mismatch version: " << version);
}

void TrialSelectPerturbed::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_(ostr);
  feasst_serialize_version(607, ostr);
}

}  // namespace feasst
