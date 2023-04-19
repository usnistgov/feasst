#include "utils/include/serialize.h"
#include "utils/include/utils.h"  // find_in_list
#include "math/include/random.h"
#include "chain/include/select_particle_pivot.h"

namespace feasst {

SelectParticlePivot::SelectParticlePivot(argtype args)
  : SelectParticlePivot(&args) {
  FEASST_CHECK_ALL_USED(args);
}

SelectParticlePivot::SelectParticlePivot(argtype * args) : TrialSelect(args) {
  class_name_ = "SelectParticlePivot";
  pivot_site_ = integer("pivot_site", args, 0);
}

class MapSelectParticlePivot {
 public:
  MapSelectParticlePivot() {
    auto obj = MakeSelectParticlePivot();
    obj->deserialize_map()["SelectParticlePivot"] = obj;
  }
};

static MapSelectParticlePivot mapper_ = MapSelectParticlePivot();

void SelectParticlePivot::precompute(System * system) {
  TrialSelect::precompute(system);
  ASSERT(is_particle_type_set(), "required particle_type as argument");
  anchor_.clear();
  anchor_.add_site(0, pivot_site_);
}

bool SelectParticlePivot::select(const Select& perturbed,
                                 System* system,
                                 Random * random) {
  const Configuration& config = system->configuration();
  const int num = config.num_particles(group_index());
  if (num <= 0) return false;
  set_probability_(1./static_cast<double>(num));
  const int index = random->uniform(0, num - 1);
  const Select& group = config.group_select(group_index());
  const int particle_index = group.particle_index(index);
  if (mobile_.num_sites() == 0) {
    mobile_ = Select(particle_index, config.select_particle(particle_index));
    mobile_.remove_sites(particle_index, {pivot_site_});
    ASSERT(mobile_.num_sites() > 0,
      "no point pivoting a particle with one site");
    ASSERT(mobile_.num_sites() == config.select_particle(particle_index).num_sites() - 1, "err");
    mobile_.resize_positions();
  } else {
    mobile_.set_particle(0, particle_index);
  }
  mobile_.load_positions(config.particles());
  anchor_.set_particle(0, particle_index);
  DEBUG("selected " << mobile_.str());
  remove_unphysical_sites(config);
  ASSERT(mobile_.num_particles() > 0, "all sites shouldn't be unphysical");
  set_mobile_original(system);
  DEBUG("selected " << mobile_.str());
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
