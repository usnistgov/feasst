#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "configuration/include/select.h"
#include "configuration/include/configuration.h"
#include "chain/include/select_branch.h"

namespace feasst {

SelectBranch::SelectBranch(argtype args) : SelectBranch(&args) {
  feasst_check_all_used(args);
}
SelectBranch::SelectBranch(argtype * args) : TrialSelectAngle(args) {
  class_name_ = "SelectBranch";
  mobile_site2_name_ = str("mobile_site2", args);
}

FEASST_MAPPER(SelectBranch,
  argtype({{"mobile_site", "2"}, {"mobile_site2", "3"},
           {"anchor_site", "1"}, {"anchor_site2", "0"}}));

void SelectBranch::precompute(System * system) {
  TrialSelectAngle::precompute(system);
  const Configuration& conf = configuration(*system);
  const int mobile_site2 = conf.site_name_to_index(mobile_site2_name_);
  get_mobile()->add_site(0, mobile_site2);
}

std::shared_ptr<TrialSelect> SelectBranch::create(std::istream& istr) const {
  return std::make_shared<SelectBranch>(istr);
}

SelectBranch::SelectBranch(std::istream& istr)
  : TrialSelectAngle(istr) {
  // ASSERT(class_name_ == "SelectBranch", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 9054 && version <= 9055, "mismatch version: " << version);
  if (version <= 9054) {
    WARN("Restart versions may be incompatible");
    int mobile_site2;
    feasst_deserialize(&mobile_site2, istr);
  }
  if (version >= 9055) {
    feasst_deserialize(&mobile_site2_name_, istr);
  }
}

void SelectBranch::serialize_select_branch_(std::ostream& ostr) const {
  serialize_trial_select_angle_(ostr);
  feasst_serialize_version(9055, ostr);
  feasst_serialize(mobile_site2_name_, ostr);
}

void SelectBranch::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_select_branch_(ostr);
}

}  // namespace feasst
