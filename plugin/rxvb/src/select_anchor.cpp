#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "configuration/include/select.h"
#include "configuration/include/particle.h"
#include "configuration/include/configuration.h"
#include "rxvb/include/select_anchor.h"

namespace feasst {

FEASST_MAPPER(SelectAnchor,);

SelectAnchor::SelectAnchor(argtype * args) : TrialSelect(args) {
  class_name_ = "SelectAnchor";
}
SelectAnchor::SelectAnchor(argtype args) : SelectAnchor(&args) {
  feasst_check_all_used(args);
}

std::shared_ptr<TrialSelect> SelectAnchor::create(std::istream& istr) const {
  return std::make_shared<SelectAnchor>(istr);
}

SelectAnchor::SelectAnchor(std::istream& istr)
  : TrialSelect(istr) {
  // ASSERT(class_name_ == "SelectAnchor", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(7556 == version, "mismatch version: " << version);
}

void SelectAnchor::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_(ostr);
  feasst_serialize_version(7556, ostr);
}

bool SelectAnchor::select(const Select& perturbed,
    System* system,
    Random * random,
    TrialSelect * previous_select) {
  if (previous_select->anchor().num_sites() == 0) return false;
  const int pindex = previous_select->anchor().particle_index(0);
  // The follow code is copied from a section of SelectParticleAVB with rxvb
  const Configuration& config = configuration(*system);
  const Particle& part = config.select_particle(pindex);
  if (mobile().num_sites() != part.num_sites()) {
    std::vector<int> sites(part.num_sites());
    for (int site = 0; site < part.num_sites(); ++site) {
      sites[site] = site;
    }
    get_mobile()->replace_indices(pindex, sites);
    get_mobile()->resize_positions();
  }
  get_mobile()->set_particle(0, pindex);
  get_mobile()->load_positions(config.particles());
  //replace_mobile(previous_select->anchor(), 0, configuration(*system));
  set_mobile_original(system);
  return true;
}

}  // namespace feasst
