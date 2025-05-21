#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/trial_select_dihedral.h"

namespace feasst {

TrialSelectDihedral::TrialSelectDihedral(argtype args) : TrialSelectDihedral(&args) {
  feasst_check_all_used(args);
}
TrialSelectDihedral::TrialSelectDihedral(argtype * args) : TrialSelectAngle(args) {
  class_name_ = "TrialSelectDihedral";
  anchor_site3_name_ = str("anchor_site3", args);
}

FEASST_MAPPER(TrialSelectDihedral,
  argtype({{"mobile_site", "0"},
      {"anchor_site", "1"}, {"anchor_site2", "2"}, {"anchor_site3", "3"}}));

void TrialSelectDihedral::precompute(System * system) {
  TrialSelectAngle::precompute(system);
  const Configuration& conf = configuration(*system);
  const int mobile_site = conf.site_name_to_index(mobile_site_name());
  const int anchor_site = conf.site_name_to_index(anchor_site_name());
  const int anchor_site2 = conf.site_name_to_index(anchor_site2_name());
  const int anchor_site3 = conf.site_name_to_index(anchor_site3_name_);
  const Particle& part = conf.particle_type(particle_type());
  const int dihedral_type = part.dihedral(mobile_site, anchor_site,
                                          anchor_site2, anchor_site3).type();
  DEBUG("mobile: " << mobile_site);
  DEBUG("anchor: " << anchor_site);
  DEBUG("anchor2: " << anchor_site2);
  DEBUG("anchor3: " << anchor_site3);
  DEBUG("dihedral_type: " << dihedral_type);
  add_or_set_property("dihedral_type", dihedral_type);
  get_anchor()->add_site(0, anchor_site3);
}

std::shared_ptr<TrialSelect> TrialSelectDihedral::create(std::istream& istr) const {
  return std::make_shared<TrialSelectDihedral>(istr);
}

TrialSelectDihedral::TrialSelectDihedral(std::istream& istr)
  : TrialSelectAngle(istr) {
  // ASSERT(class_name_ == "TrialSelectDihedral", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 8954 && version <= 8955, "mismatch version: " << version);
  if (version <= 8954) {
    WARN("Restart versions may be incompatible");
    int anchor_site3;
    feasst_deserialize(&anchor_site3, istr);

  }
  if (version >= 8955) {
    feasst_deserialize(&anchor_site3_name_, istr);
  }
}

void TrialSelectDihedral::serialize_trial_select_dihedral_(std::ostream& ostr) const {
  serialize_trial_select_angle_(ostr);
  feasst_serialize_version(8955, ostr);
  feasst_serialize(anchor_site3_name_, ostr);
}

void TrialSelectDihedral::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_dihedral_(ostr);
}

}  // namespace feasst
