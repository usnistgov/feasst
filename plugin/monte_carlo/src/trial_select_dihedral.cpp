#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_dihedral.h"

namespace feasst {

TrialSelectDihedral::TrialSelectDihedral(argtype args) : TrialSelectDihedral(&args) {
  FEASST_CHECK_ALL_USED(args);
}
TrialSelectDihedral::TrialSelectDihedral(argtype * args) : TrialSelectAngle(args) {
  class_name_ = "TrialSelectDihedral";
  anchor_site3_ = integer("anchor_site3", args);
}

class MapTrialSelectDihedral {
 public:
  MapTrialSelectDihedral() {
    auto obj = MakeTrialSelectDihedral({{"mobile_site", "0"},
      {"anchor_site", "1"}, {"anchor_site2", "2"}, {"anchor_site3", "3"}});
    obj->deserialize_map()["TrialSelectDihedral"] = obj;
  }
};

static MapTrialSelectDihedral mapper_ = MapTrialSelectDihedral();

void TrialSelectDihedral::precompute(System * system) {
  TrialSelectAngle::precompute(system);
  const Particle& part = configuration(*system).particle_type(particle_type());
  const int dihedral_type = part.dihedral(mobile_site(), anchor_site(),
                                          anchor_site2(), anchor_site3_).type();
  DEBUG("mobile: " << mobile_site());
  DEBUG("anchor: " << anchor_site());
  DEBUG("anchor2: " << anchor_site2());
  DEBUG("anchor3: " << anchor_site3_);
  DEBUG("dihedral_type: " << dihedral_type);
  add_or_set_property("dihedral_type", dihedral_type);
  anchor_.add_site(0, anchor_site3_);
}

std::shared_ptr<TrialSelect> TrialSelectDihedral::create(std::istream& istr) const {
  return std::make_shared<TrialSelectDihedral>(istr);
}

TrialSelectDihedral::TrialSelectDihedral(std::istream& istr)
  : TrialSelectAngle(istr) {
  // ASSERT(class_name_ == "TrialSelectDihedral", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(8954 == version, "mismatch version: " << version);
  feasst_deserialize(&anchor_site3_, istr);
}

void TrialSelectDihedral::serialize_trial_select_dihedral_(std::ostream& ostr) const {
  serialize_trial_select_angle_(ostr);
  feasst_serialize_version(8954, ostr);
  feasst_serialize(anchor_site3_, ostr);
}

void TrialSelectDihedral::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_dihedral_(ostr);
}

}  // namespace feasst
