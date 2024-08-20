#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
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
  feasst_check_all_used(args);
}
SelectReptate::SelectReptate(argtype * args) : SelectEndSegment(args) {
  ASSERT(max_length() == 1, "requires max_length(" << max_length() << ") of 1");
  class_name_ = "SelectReptate";
}

void SelectReptate::precompute(System * system) {
  SelectEndSegment::precompute(system);
  { // select bond
    const Particle& part = configuration(*system).particle_types().particle(particle_type());
    const int bond_type = part.bond(0, 1).type();  // HWH assume constant bond length
    add_or_set_property("bond_type", bond_type);
  }
  get_anchor()->clear();
  get_anchor()->add_site(0, 0);
  bonded_to_.clear();
  bonded_to_.add_site(0, 0);
}

void SelectReptate::update_anchor(const bool is_endpoint_beginning,
  const System * system) {
  const int particle_index = mobile().particle_indices()[0];
  const Configuration& config = configuration(*system);
  const Particle& particle = config.select_particle(particle_index);
  int anchor_index = 0;
  int site_bonded_to = particle.num_sites() - 2;
  DEBUG("is_endpoint_beginning " << is_endpoint_beginning);
  if (is_endpoint_beginning) {
    anchor_index = particle.num_sites() - 1;
    site_bonded_to = 1;
  }
  // for the old configuration, set the anchor to the old bond.
  get_anchor()->set_site(0, 0, anchor_index);
  get_anchor()->set_particle(0, particle_index);
  ASSERT(bonded_to_.replace_indices(particle_index, {site_bonded_to}),
    "bonded_to_ wasn't initialized to proper size on precompute");
}

void SelectReptate::mid_stage() {
  // exclude the anchor from interactions.
  // include interactions with site that use to be bonded
  get_mobile()->set_new_bond(anchor());
  get_mobile()->set_old_bond(bonded_to_);
}

}  // namespace feasst
