#include "chain/include/trial_protonation.h"

namespace feasst {

class MapTrialProtonation {
 public:
  MapTrialProtonation() {
    auto obj = MakeTrialProtonation({
      {"reactant_type", "0"},
      {"reactant_site_type", "0"},
      {"new_site_type", "1"},
      {"remove_type", "1"}
    });
    obj->deserialize_map()["TrialProtonation"] = obj;
  }
};

static MapTrialProtonation mapper_ = MapTrialProtonation();

std::shared_ptr<Trial> TrialProtonation::create(std::istream& istr) const {
  return std::make_shared<TrialProtonation>(istr);
}

TrialProtonation::TrialProtonation(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialProtonation", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(179 == version, "mismatch version: " << version);
}

void TrialProtonation::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(179, ostr);
}

}  // namespace feasst
