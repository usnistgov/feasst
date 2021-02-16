#include "utils/include/serialize.h"
#include "chain/include/select_branch.h"

namespace feasst {

SelectBranch::SelectBranch(argtype args) : SelectBranch(&args) {
  check_all_used(args);
}
SelectBranch::SelectBranch(argtype * args) : TrialSelectAngle(args) {
  class_name_ = "SelectBranch";
  mobile_site2_ = integer("mobile_site2", args);
}

class MapSelectBranch {
 public:
  MapSelectBranch() {
    auto obj = MakeSelectBranch({{"mobile_site", "2"}, {"mobile_site2", "3"},
                                 {"anchor_site", "1"}, {"anchor_site2", "0"}});
    obj->deserialize_map()["SelectBranch"] = obj;
  }
};

static MapSelectBranch mapper_ = MapSelectBranch();

void SelectBranch::precompute(System * system) {
  TrialSelectAngle::precompute(system);
  mobile_.add_site(0, mobile_site2_);
}

std::shared_ptr<TrialSelect> SelectBranch::create(std::istream& istr) const {
  return std::make_shared<SelectBranch>(istr);
}

SelectBranch::SelectBranch(std::istream& istr)
  : TrialSelectAngle(istr) {
  // ASSERT(class_name_ == "SelectBranch", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(9054 == version, "mismatch version: " << version);
  feasst_deserialize(&mobile_site2_, istr);
}

void SelectBranch::serialize_select_branch_(std::ostream& ostr) const {
  serialize_trial_select_angle_(ostr);
  feasst_serialize_version(9054, ostr);
  feasst_serialize(mobile_site2_, ostr);
}

void SelectBranch::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_select_branch_(ostr);
}

}  // namespace feasst
