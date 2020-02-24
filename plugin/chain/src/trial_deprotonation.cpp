#include "chain/include/trial_deprotonation.h"

namespace feasst {

class MapTrialDeprotonation {
 public:
  MapTrialDeprotonation() {
    auto obj = MakeTrialDeprotonation({
      {"reactant_type", "0"},
      {"reactant_site_type", "0"},
      {"new_site_type", "1"},
      {"add_type", "1"}
    });
    obj->deserialize_map()["TrialDeprotonation"] = obj;
  }
};

static MapTrialDeprotonation mapper_ = MapTrialDeprotonation();

std::shared_ptr<Trial> TrialDeprotonation::create(std::istream& istr) const {
  return std::make_shared<TrialDeprotonation>(istr);
}

TrialDeprotonation::TrialDeprotonation(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialDeprotonation", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(368 == version, "mismatch version: " << version);
}

void TrialDeprotonation::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(368, ostr);
}

}  // namespace feasst
