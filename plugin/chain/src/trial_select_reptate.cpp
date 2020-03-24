#include "utils/include/serialize.h"
#include "chain/include/trial_select_reptate.h"

namespace feasst {

class MapTrialSelectReptate {
 public:
  MapTrialSelectReptate() {
    auto obj = MakeTrialSelectReptate({{"max_length", "1"}});
    obj->deserialize_map()["TrialSelectReptate"] = obj;
  }
};

static MapTrialSelectReptate mapper_ = MapTrialSelectReptate();

std::shared_ptr<TrialSelect> TrialSelectReptate::create(std::istream& istr) const {
  return std::make_shared<TrialSelectReptate>(istr);
}

TrialSelectReptate::TrialSelectReptate(std::istream& istr)
  : TrialSelectEndSegment(istr) {
  // ASSERT(class_name_ == "TrialSelectReptate", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(812 == version, "mismatch version: " << version);
  feasst_deserialize_fstobj(&bonded_to_, istr);
}

void TrialSelectReptate::serialize_trial_select_reptate_(std::ostream& ostr) const {
  serialize_trial_select_end_segment_(ostr);
  feasst_serialize_version(812, ostr);
  feasst_serialize_fstobj(bonded_to_, ostr);
}

void TrialSelectReptate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_reptate_(ostr);
}

TrialSelectReptate::TrialSelectReptate(const argtype& args)
  : TrialSelectEndSegment(args) {
  ASSERT(max_length() == 1,
    "requires max_length(" << max_length() << ") of 1");
  class_name_ = "TrialSelectReptate";
}

void TrialSelectReptate::precompute(System * system) {
  TrialSelectEndSegment::precompute(system);
  anchor_.clear();
  anchor_.add_site(0, 0);
  bonded_to_.clear();
  bonded_to_.add_site(0, 0);
}

void TrialSelectReptate::update_anchor(const bool is_endpoint_beginning,
  const System * system) {
  const int particle_index = mobile_.particle_indices()[0];
  const Configuration& config = system->configuration();
  const Particle& particle = config.select_particle(particle_index);
  int anchor_index = 0;
  int site_bonded_to = particle.num_sites() - 2;
  DEBUG("is_endpoint_beginning " << is_endpoint_beginning);
  if (is_endpoint_beginning) {
    anchor_index = particle.num_sites() - 1;
    site_bonded_to = 1;
  }
  // for the old configuration, set the anchor to the old bond.
  anchor_.set_site(0, 0, anchor_index);
  anchor_.set_particle(0, particle_index);
  ASSERT(bonded_to_.replace_indices(particle_index, {site_bonded_to}),
    "bonded_to_ wasn't initialized to proper size on precompute");
}

}  // namespace feasst
