#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "utils/include/utils.h"  // find_in_list
#include "math/include/random.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "chain/include/select_particle_pivot.h"

namespace feasst {

SelectParticlePivot::SelectParticlePivot(argtype args)
  : SelectParticlePivot(&args) {
  feasst_check_all_used(args);
}

SelectParticlePivot::SelectParticlePivot(argtype * args) : TrialSelect(args) {
  class_name_ = "SelectParticlePivot";
  pivot_site_ = integer("pivot_site", args, 0);
}

FEASST_MAPPER(SelectParticlePivot,);

void SelectParticlePivot::precompute(System * system) {
  TrialSelect::precompute(system);
  ASSERT(is_particle_type_set(), "required particle_type as argument");
  ASSERT(particle_type() >= 0, "particle_type required and must be >= 0");
  get_anchor()->clear();
  get_anchor()->add_site(0, pivot_site_);
}

bool SelectParticlePivot::select(const Select& perturbed,
                                 System* system,
                                 Random * random) {
  const Configuration& config = configuration(*system);
  const int num = config.num_particles(group_index());
  if (num <= 0) return false;
  set_probability_(1./static_cast<double>(num));
  const int index = random->uniform(0, num - 1);
  const Select& group = config.group_select(group_index());
  const int particle_index = group.particle_index(index);
  if (mobile().num_sites() == 0) {
    set_mobile(Select(particle_index, config.select_particle(particle_index)));
    get_mobile()->remove_sites(particle_index, {pivot_site_});
    ASSERT(mobile().num_sites() > 0,
      "no point pivoting a particle with one site");
    ASSERT(mobile().num_sites() == config.select_particle(particle_index).num_sites() - 1, "err");
    get_mobile()->resize_positions();
  } else {
    get_mobile()->set_particle(0, particle_index);
  }
  get_mobile()->load_positions(config.particles());
  get_anchor()->set_particle(0, particle_index);
  DEBUG("selected " << mobile().str());
  remove_unphysical_sites(config);
  ASSERT(mobile().num_particles() > 0, "all sites shouldn't be unphysical");
  set_mobile_original(system);
  DEBUG("selected " << mobile().str());
  return true;
}

std::shared_ptr<TrialSelect> SelectParticlePivot::create(std::istream& istr) const {
  return std::make_shared<SelectParticlePivot>(istr);
}

SelectParticlePivot::SelectParticlePivot(std::istream& istr)
  : TrialSelect(istr) {
  // ASSERT(class_name_ == "SelectParticlePivot", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(6841 == version, "mismatch version: " << version);
  feasst_deserialize(&pivot_site_, istr);
}

void SelectParticlePivot::serialize_select_particle_pivot_(std::ostream& ostr) const {
  serialize_trial_select_(ostr);
  feasst_serialize_version(6841, ostr);
  feasst_serialize(pivot_site_, ostr);
}

void SelectParticlePivot::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_select_particle_pivot_(ostr);
}

}  // namespace feasst
