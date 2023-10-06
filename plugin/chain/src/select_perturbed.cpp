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

SelectPerturbed::SelectPerturbed(argtype * args) : TrialSelect(args) {
  class_name_ = "SelectPerturbed";
}
SelectPerturbed::SelectPerturbed(argtype args) : SelectPerturbed(&args) {
  FEASST_CHECK_ALL_USED(args);
}

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

bool SelectPerturbed::select(const Select& perturbed, System* system, Random * random) {
  if (perturbed.num_sites() == 0) return false;
  replace_mobile(perturbed, 0, configuration(*system));
  set_mobile_original(system);
  return true;
}

}  // namespace feasst
