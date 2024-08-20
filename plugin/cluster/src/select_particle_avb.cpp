#include "utils/include/serialize.h"
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/utils.h"
#include "math/include/utils_math.h"
#include "math/include/random.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/neighbor_criteria.h"
#include "configuration/include/configuration.h"
#include "configuration/include/domain.h"
#include "system/include/system.h"
#include "system/include/energy_map.h"
#include "cluster/include/select_particle_avb.h"

namespace feasst {

SelectParticleAVB::SelectParticleAVB(argtype args) : SelectParticleAVB(&args) {
  feasst_check_all_used(args);
}
SelectParticleAVB::SelectParticleAVB(argtype * args) : TrialSelect(args) {
  class_name_ = "SelectParticleAVB";
  neighbor_ = integer("neighbor_index", args, 0);
  site_index_ = integer("site", args, 0);
  get_mobile()->clear();
  get_mobile()->add_site(0, site_index_);
  grand_canonical_ = boolean("grand_canonical", args);
  inside_ = boolean("inside", args, true);
  is_second_target_ = boolean("second_target", args, false);

  // initialize select_target_
  argtype target_args;
  target_args.insert({"load_coordinates", "true"});
  target_args.insert({"particle_type",
                      str("target_particle_type", args, "0")});
  target_args.insert({"site", str("target_site", args, "0")});
  select_target_ = TrialSelectParticle(target_args);

  // initialize select_mobile_
  argtype mobile_args;
  mobile_args.insert({"load_coordinates", "true"});
  if (is_particle_type_set()) {
    mobile_args.insert({"particle_type", str(particle_type())});
  }
  mobile_args.insert({"site", str(site_index_)});
  select_mobile_ = TrialSelectParticle(mobile_args);

  ASSERT(!used("group_index", *args), "group not implemented with AVB");
}

class MapSelectParticleAVB {
 public:
  MapSelectParticleAVB() {
    auto obj = MakeSelectParticleAVB({{"grand_canonical", "true"}});
    obj->deserialize_map()["SelectParticleAVB"] = obj;
  }
};

static MapSelectParticleAVB mapper_ = MapSelectParticleAVB();

void SelectParticleAVB::precompute(System * system) {
  TrialSelect::precompute(system);
  select_target_.precompute(system);
  select_mobile_.precompute(system);
  get_anchor()->clear();
  get_anchor()->add_site(0, select_target_.site());
}

bool SelectParticleAVB::select(const Select& perturbed,
                               System * system,
                               Random * random) {
  const Configuration& config = system->configuration();
  if ( (is_ghost() && (config.num_particles() < 1)) ||
       (!is_ghost() && (config.num_particles() < 2)) ) {
    DEBUG("avb not possible if n=" << config.num_particles() <<
      " and ghost:" << is_ghost());
    return false;
  }
  DEBUG("target sel type: " << select_target_.particle_type());
  const int num = select_target_.random_particle(config, &target_, random);
  DEBUG("num " << num);
  if (num <= 0) {
    DEBUG("avb not possible");
    return false;
  }
  DEBUG("target: " << target_.str());
  DEBUG("target type: " << config.select_particle(target_.particle_index(0)).type());

  // set anchor
  if (is_second_target_) {
    const int num2 = select_target_.random_particle(
      config, &target_, &second_target_, random);
    if (num2 <= 0) return false;
    get_anchor()->set_particle(0, second_target_.particle_index(0));
    DEBUG("second_target " << second_target_.str());
  } else {
    get_anchor()->set_particle(0, target_.particle_index(0));
  }
  ASSERT(target_.num_sites() == 1, "Error");
  // HWH update with configuration_index_
  const NeighborCriteria& neighbor = system->neighbor_criteria(neighbor_, 0);
  map_(*system, neighbor_).neighbors(
    neighbor,
    config,
    target_.particle_index(0),
    target_.site_index(0, 0),
    site_index_,
    &neighbors_);
  const int num_neighbors = static_cast<int>(neighbors_.num_sites());
  DEBUG("neighbors: " << neighbors_.str());

  // Initialize mobile
  int num_out = 0.;
  if (!is_ghost() && inside_) {
    DEBUG("inside. num_neighbors: " << num_neighbors);
    if (num_neighbors <= 0) return false;
    if (is_second_target_) {
      if (second_target_.is_overlap(neighbors_)) {
        return false;
      }
    }
    get_mobile()->set_particle(0,
      random->const_element(neighbors_.particle_indices()));
    DEBUG("mobile " << mobile().str());
  } else if (!inside_) {
    // only select that is outside: AVB2 out->in
    DEBUG("outside");
    neighbors_.add(target_); // add target to neighbors to exclude.
    num_out = select_mobile_.random_particle(config,
      &neighbors_,
      get_mobile(),
      random);
    DEBUG("num_out: " << num_out);
    if (num_out == 0) return false;
  }
  DEBUG("loading coordinates " << mobile().num_particles());
  get_mobile()->load_positions(config.particles());

  // precompute volume terms
  const double volume = config.domain().volume();
  const double volume_av = neighbor.volume(config.dimension());
  const double volume_out = volume - volume_av;
  ASSERT(volume_out > 0, "AV volume: " << volume_av << " is too large "
    << "for total volume: " << volume);

  bool target_mobile_same_type_ = true;
  if (is_particle_type_set()) {
    if (select_target_.particle_type() != particle_type()) {
      target_mobile_same_type_ = false;
    }
  }

  // Assign the probabilities for acceptance criteria.
  // explicitly consider all valid combinations the following booleans:
  // is_ghost, grand_canoical, inside and second_target.
  // There are five valid combinations which represent GCE add, GCE rm,
  // AVB2 in->out, AVB2 out->in and AVB4 in->in.

  // GCE add
  if (is_ghost() && grand_canonical_ && inside_ && !is_second_target_) {
    select_mobile_.ghost_particle(
      system->get_configuration(), &empty_, get_mobile());
    DEBUG("num_neighbors " << num_neighbors);
    set_probability_(volume_av/static_cast<double>(num_neighbors + 1));
    DEBUG("target_mobile_same_type_ " << target_mobile_same_type_);
    if (target_mobile_same_type_) {
      DEBUG("prob before: " << probability());
      DEBUG("num " << num);
      // ghost will be added during perturb. Not yet added, so n/(n+1) factor.
      set_probability_(probability()
        *static_cast<double>(num)/static_cast<double>(num + 1));
    }

  // GCE remove
  } else if (!is_ghost() && grand_canonical_ && inside_ && !is_second_target_) {
    set_probability_(static_cast<double>(num_neighbors)/volume_av);
    if (target_mobile_same_type_) {
      set_probability_(probability()
        *static_cast<double>(num)/static_cast<double>(num - 1));
    }

  // AVB2 in->out
  } else if (!is_ghost() && !grand_canonical_ && inside_ && !is_second_target_) {
    //DEBUG("AVB2 in->out");
    // compute num_out
    ASSERT(num_out == 0, "num_out from above should be zero");
    int num_tot_tmp;
    if (is_particle_type_set()) {
      num_tot_tmp = config.num_particles_of_type(particle_type());
    } else {
      num_tot_tmp = config.num_particles();
    }
    num_out = num_tot_tmp - num_neighbors;
    if (target_mobile_same_type_) {
      --num_out;
    }
    // ComputeAVB2 will add p_bias term
    set_probability_(static_cast<double>(num_neighbors)
      /static_cast<double>(num_out + 1)
      *volume_out/volume_av);

  // AVB2 out->in
  } else if (!is_ghost() && !grand_canonical_ && !inside_ && !is_second_target_) {
    // ComputeAVB2 will add p_bias term
    set_probability_(static_cast<double>(num_out)
      /static_cast<double>(num_neighbors + 1)
      *volume_av/volume_out);

  // AVB4 in->in
  } else if (!is_ghost() && !grand_canonical_ && inside_ && is_second_target_) {
    // obtain the number of neighbors in the second target
    map_(*system, neighbor_).neighbors(
      neighbor,
      config,
      second_target_.particle_index(0),
      second_target_.site_index(0, 0),
      site_index_,
      &neighbors_);
    DEBUG("avb4 neighbors: " << neighbors_.str());
    DEBUG("avb4 neighbors in second target: " << neighbors_.num_particles());
    int mobile_not_in_both_targets_ = 1;
    DEBUG("mobile part ind " << mobile().particle_index(0));
    DEBUG("mobile k neighs: " << neighbors_.str());
    DEBUG("mobile: " << mobile().str());
    DEBUG("second_target: " << second_target_.str());
    if (find_in_list(mobile().particle_index(0),
                     neighbors_.particle_indices())) {
      mobile_not_in_both_targets_ = 0;
    }
    DEBUG("mobile_not_in_both_targets_ " << mobile_not_in_both_targets_);
    set_probability_(static_cast<double>(num_neighbors)
      /static_cast<double>(neighbors_.num_particles() +
                           mobile_not_in_both_targets_));

  // invalid
  } else {
    FATAL("unrecognized combination of the is_ghost: " << is_ghost()
      << " grand_canonical: " << grand_canonical_ << " inside: " << inside_
      << " second_target: " << is_second_target_);
  }

  ASSERT(mobile().particle_index(0) != target_.particle_index(0),
    "mobile particle: " << mobile().particle_index(0) <<
    " should not be equal to target: " << target_.particle_index(0));
  DEBUG("probability: " << probability());
  DEBUG("mobile: " << mobile().str());
//  printable_["cluster_size"].accumulate(mobile().num_particles());
  remove_unphysical_sites(config);
  set_mobile_original(system);
  return true;
}

std::shared_ptr<TrialSelect> SelectParticleAVB::create(std::istream& istr) const {
  return std::make_shared<SelectParticleAVB>(istr);
}

SelectParticleAVB::SelectParticleAVB(std::istream& istr)
  : TrialSelect(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(239 == version, "mismatch version: " << version);
  feasst_deserialize(&neighbor_, istr);
  feasst_deserialize(&site_index_, istr);
  feasst_deserialize(&grand_canonical_, istr);
  feasst_deserialize(&inside_, istr);
  feasst_deserialize(&is_second_target_, istr);
  feasst_deserialize_fstobj(&select_target_, istr);
  feasst_deserialize_fstobj(&select_mobile_, istr);
}

void SelectParticleAVB::serialize_select_particle_avb_(
    std::ostream& ostr) const {
  serialize_trial_select_(ostr);
  feasst_serialize_version(239, ostr);
  feasst_serialize(neighbor_, ostr);
  feasst_serialize(site_index_, ostr);
  feasst_serialize(grand_canonical_, ostr);
  feasst_serialize(inside_, ostr);
  feasst_serialize(is_second_target_, ostr);
  feasst_serialize_fstobj(select_target_, ostr);
  feasst_serialize_fstobj(select_mobile_, ostr);
}

void SelectParticleAVB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_select_particle_avb_(ostr);
}

}  // namespace feasst
