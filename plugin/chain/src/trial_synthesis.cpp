#include "chain/include/trial_synthesis.h"

namespace feasst {

class MapTrialSynthesis {
 public:
  MapTrialSynthesis() {
    auto obj = MakeTrialSynthesis({
      {"reactant_type", "0"},
      {"reactant_site_type", "0"},
      {"new_site_type", "1"},
      {"product_type", "1"}
    });
    obj->deserialize_map()["TrialSynthesis"] = obj;
  }
};

static MapTrialSynthesis mapper_ = MapTrialSynthesis();

std::shared_ptr<Trial> TrialSynthesis::create(std::istream& istr) const {
  return std::make_shared<TrialSynthesis>(istr);
}

TrialSynthesis::TrialSynthesis(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialSynthesis", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(368 == version, "mismatch version: " << version);
}

void TrialSynthesis::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(368, ostr);
}

}  // namespace feasst
