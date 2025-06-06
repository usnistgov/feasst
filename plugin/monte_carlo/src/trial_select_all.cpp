#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "utils/include/utils.h"  // find_in_list
#include "math/include/random.h"
#include "configuration/include/select.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/trial_select_all.h"

namespace feasst {

TrialSelectAll::TrialSelectAll(argtype * args) : TrialSelect(args) {
  class_name_ = "TrialSelectAll";
}
TrialSelectAll::TrialSelectAll(argtype args) : TrialSelectAll(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(TrialSelectAll,);

bool TrialSelectAll::select(const Select& perturbed,
                                  System* system,
                                  Random * random,
                                  TrialSelect * previous_select) {
  set_probability_(1.);
  const Configuration& config = configuration(*system);
  set_mobile(config.selection_of_all());
  DEBUG("selected " << mobile().str());
  remove_unphysical_sites(config);
  set_mobile_original(system);
  DEBUG("selected " << mobile().str());
  return true;
}

std::shared_ptr<TrialSelect> TrialSelectAll::create(std::istream& istr) const {
  return std::make_shared<TrialSelectAll>(istr);
}

TrialSelectAll::TrialSelectAll(std::istream& istr)
  : TrialSelect(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(1759 == version, "mismatch version: " << version);
}

void TrialSelectAll::serialize_trial_select_all_(std::ostream& ostr) const {
  serialize_trial_select_(ostr);
  feasst_serialize_version(1759, ostr);
}

void TrialSelectAll::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_all_(ostr);
}

}  // namespace feasst
