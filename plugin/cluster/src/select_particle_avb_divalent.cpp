#include "utils/include/serialize.h"
#include "utils/include/utils.h"
#include "math/include/utils_math.h"
#include "math/include/random.h"
#include "configuration/include/domain.h"
#include "cluster/include/select_particle_avb_divalent.h"

namespace feasst {

SelectParticleAVBDivalent::SelectParticleAVBDivalent(argtype args)
  : TrialSelect(&args) {
  class_name_ = "SelectParticleAVBDivalent";
  neighbor_ = integer("neighbor_index", &args, 0);

  // initialize select_mobile_
  argtype mobile_args;
  mobile_args.insert({"load_coordinates", "true"});
  mobile_args.insert({"particle_type", str(particle_type())});
  mobile_args.insert({"site", str("site_index", &args, "0")});
  mobile_args.insert({"ghost", str("ghost", &args)});
  select_mobile_ = TrialSelectParticle(mobile_args);
  DEBUG(select_mobile_.particle_type());
  mobile_.clear();
  mobile_.add_site(0, select_mobile_.site());

  anchor_.clear();
  anchor_.add_site(0,
    integer("target_site_index", &args, 0));
  ASSERT(!used("group_index", args), "group not implemented with AVB");
  check_all_used(args);
}

class MapSelectParticleAVBDivalent {
 public:
  MapSelectParticleAVBDivalent() {
    auto obj = MakeSelectParticleAVBDivalent(
      {{"ghost", "true"}, {"particle_type", "0"}});
    obj->deserialize_map()["SelectParticleAVBDivalent"] = obj;
  }
};

static MapSelectParticleAVBDivalent mapper_ = MapSelectParticleAVBDivalent();

bool SelectParticleAVBDivalent::select(const Select& perturbed,
                               System * system,
                               Random * random) {
  const Configuration& config = system->configuration();
  ASSERT(perturbed.num_particles() >= 1, "first stage should have completed.");
  anchor_.set_particle(0, perturbed.particle_index(0));
  if (select_mobile_.is_ghost()) {
    select_mobile_.ghost_particle(
      system->get_configuration(),
      const_cast<const Select*>(&perturbed),
      &mobile_);
  } else {
    const NeighborCriteria& neighbor = system->neighbor_criteria(neighbor_);
    map_(*system, neighbor_).neighbors(
      neighbor,
      config,
      anchor_.particle_index(0),
      anchor_.site_index(0, 0),
      mobile_.site_index(0, 0),
      &neighbors_);
    neighbors_.remove(perturbed);
    const int num_neigh = static_cast<int>(neighbors_.num_sites());
    if (num_neigh == 0) return false;
    DEBUG("num neigh " << num_neigh << " : " << neighbors_.str());
    mobile_.set_particle(0,
      random->const_element(neighbors_.particle_indices()));
    const double volume_av = neighbor.volume(config.dimension());
    set_probability_(num_neigh/volume_av);
  }
  DEBUG("probability: " << probability());
  DEBUG("mobile: " << mobile_.str());
  remove_unphysical_sites(config);
  mobile_original_ = mobile_;
  return true;
}

std::shared_ptr<TrialSelect> SelectParticleAVBDivalent::create(std::istream& istr) const {
  return std::make_shared<SelectParticleAVBDivalent>(istr);
}

SelectParticleAVBDivalent::SelectParticleAVBDivalent(std::istream& istr)
  : TrialSelect(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(2385 == version, "mismatch version: " << version);
  feasst_deserialize(&neighbor_, istr);
  feasst_deserialize_fstobj(&select_mobile_, istr);
}

void SelectParticleAVBDivalent::serialize_select_particle_avb_divalent_(
    std::ostream& ostr) const {
  serialize_trial_select_(ostr);
  feasst_serialize_version(2385, ostr);
  feasst_serialize(neighbor_, ostr);
  feasst_serialize_fstobj(select_mobile_, ostr);
}

void SelectParticleAVBDivalent::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_select_particle_avb_divalent_(ostr);
}

}  // namespace feasst
