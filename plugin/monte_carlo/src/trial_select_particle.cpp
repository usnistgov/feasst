#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

class MapTrialSelectParticle {
 public:
  MapTrialSelectParticle() {
    auto obj = MakeTrialSelectParticle();
    obj->deserialize_map()["TrialSelectParticle"] = obj;
  }
};

static MapTrialSelectParticle mapper_ = MapTrialSelectParticle();

std::shared_ptr<TrialSelect> TrialSelectParticle::create(std::istream& istr) const {
  return std::make_shared<TrialSelectParticle>(istr);
}

TrialSelectParticle::TrialSelectParticle(std::istream& istr)
  : TrialSelect(istr) {
  // ASSERT(class_name_ == "TrialSelectParticle", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(760 == version, "mismatch version: " << version);
  feasst_deserialize(&load_coordinates_, istr);
  feasst_deserialize(&site_, istr);
  feasst_deserialize(&site_vec_, istr);
}

void TrialSelectParticle::serialize_trial_select_particle_(std::ostream& ostr) const {
  serialize_trial_select_(ostr);
  feasst_serialize_version(760, ostr);
  feasst_serialize(load_coordinates_, ostr);
  feasst_serialize(site_, ostr);
  feasst_serialize(site_vec_, ostr);
}

void TrialSelectParticle::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_particle_(ostr);
}

}  // namespace feasst
