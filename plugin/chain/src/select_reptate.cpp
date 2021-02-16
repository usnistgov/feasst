#include "utils/include/serialize.h"
#include "chain/include/select_reptate.h"

namespace feasst {

class MapSelectReptate {
 public:
  MapSelectReptate() {
    auto obj = MakeSelectReptate({{"max_length", "1"}});
    obj->deserialize_map()["SelectReptate"] = obj;
  }
};

static MapSelectReptate mapper_ = MapSelectReptate();

std::shared_ptr<TrialSelect> SelectReptate::create(std::istream& istr) const {
  return std::make_shared<SelectReptate>(istr);
}

SelectReptate::SelectReptate(std::istream& istr)
  : SelectEndSegment(istr) {
  // ASSERT(class_name_ == "SelectReptate", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(812 == version, "mismatch version: " << version);
  feasst_deserialize_fstobj(&bonded_to_, istr);
}

void SelectReptate::serialize_select_reptate_(std::ostream& ostr) const {
  serialize_select_end_segment_(ostr);
  feasst_serialize_version(812, ostr);
  feasst_serialize_fstobj(bonded_to_, ostr);
}

void SelectReptate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_select_reptate_(ostr);
}

SelectReptate::SelectReptate(argtype args) : SelectReptate(&args) {
  check_all_used(args);
}
SelectReptate::SelectReptate(argtype * args) : SelectEndSegment(args) {
  ASSERT(max_length() == 1, "requires max_length(" << max_length() << ") of 1");
  class_name_ = "SelectReptate";
}

void SelectReptate::precompute(System * system) {
  SelectEndSegment::precompute(system);
  anchor_.clear();
  anchor_.add_site(0, 0);
  bonded_to_.clear();
  bonded_to_.add_site(0, 0);
}

void SelectReptate::update_anchor(const bool is_endpoint_beginning,
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
