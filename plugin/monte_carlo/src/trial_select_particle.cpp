#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

TrialSelectParticle::TrialSelectParticle(const argtype& args)
  : TrialSelect(args) {
  class_name_ = "TrialSelectParticle";
  Arguments args_(args);
  args_.dont_check();
  load_coordinates_ = args_.key("load_coordinates").dflt("true").boolean();

  // parse site
  site_ = args_.key("site").dflt("-1").integer();
  if (site_ != -1) {
    mobile_.clear();
    mobile_.add_site(0, site_);
    site_vec_ =  {site_};
  }

  set_ghost(args_.key("ghost").dflt("false").boolean());
}

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
