#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_angle.h"

namespace feasst {

TrialSelectAngle::TrialSelectAngle(const argtype& args)
  : TrialSelectBond(args) {
  class_name_ = "TrialSelectAngle";
  Arguments args_(args);
  args_.dont_check();
  anchor_site2_ = args_.key("anchor_site2").integer();
}

class MapTrialSelectAngle {
 public:
  MapTrialSelectAngle() {
    auto obj = MakeTrialSelectAngle({{"mobile_site", "2"}, {"anchor_site", "1"}, {"anchor_site2", "0"}});
    obj->deserialize_map()["TrialSelectAngle"] = obj;
  }
};

static MapTrialSelectAngle mapper_ = MapTrialSelectAngle();

void TrialSelectAngle::precompute(System * system) {
  TrialSelectBond::precompute(system);
  const Particle& part = system->configuration().particle_type(particle_type());
  const int angle_type = part.angle(mobile_site(),
                                    anchor_site(),
                                    anchor_site2_).type();
  DEBUG("mobile: " << mobile_site());
  DEBUG("anchor: " << anchor_site());
  DEBUG("anchor2: " << anchor_site2_);
  DEBUG("angle_type: " << angle_type);
  add_or_set_property("angle_type", angle_type);
  anchor_.add_site(0, anchor_site2_);
}

std::shared_ptr<TrialSelect> TrialSelectAngle::create(std::istream& istr) const {
  return std::make_shared<TrialSelectAngle>(istr);
}

TrialSelectAngle::TrialSelectAngle(std::istream& istr)
  : TrialSelectBond(istr) {
  // ASSERT(class_name_ == "TrialSelectAngle", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(666 == version, "mismatch version: " << version);
  feasst_deserialize(&anchor_site2_, istr);
}

void TrialSelectAngle::serialize_trial_select_angle_(std::ostream& ostr) const {
  serialize_trial_select_bond_(ostr);
  feasst_serialize_version(666, ostr);
  feasst_serialize(anchor_site2_, ostr);
}

void TrialSelectAngle::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_angle_(ostr);
}

}  // namespace feasst
