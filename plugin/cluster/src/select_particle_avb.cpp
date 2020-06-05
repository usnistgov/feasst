#include "utils/include/serialize.h"
#include "utils/include/utils.h"
#include "math/include/utils_math.h"
#include "math/include/random.h"
#include "configuration/include/domain.h"
#include "cluster/include/select_particle_avb.h"

namespace feasst {

SelectParticleAVB::SelectParticleAVB(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args)
  : TrialSelect(args) {
  class_name_ = "SelectParticleAVB";
  neighbor_criteria_ = neighbor_criteria;
  site_index_ = args_.key("site_index").dflt("0").integer();
  mobile_.clear();
  mobile_.add_site(0, 0);
  grand_canonical_ = args_.key("grand_canonical").boolean();
  inside_ = args_.key("inside").dflt("true").boolean();
  is_second_target_ = args_.key("second_target").dflt("false").boolean();

// HWH ghost isn't set by PerturbAdd::precompute
//  set_ghost(args_.key("ghost").dflt("false").boolean());
//  if (is_ghost()) ASSERT(grand_canonical_, "ghost requires grand_canonical");
//  DEBUG("ghost? " << is_ghost());

  // initialize select_target_
  argtype target_args;
//  if (is_ghost()) {
    target_args.insert({"load_coordinates", "true"});
//  } else {
//    target_args.insert({"load_coordinates", "false"});
//  }
  target_args.insert({"particle_type",
                      args_.key("target_particle_type").dflt("0").str()});
  target_args.insert({"site", args_.key("target_site_index").dflt("0").str()});
  select_target_ = TrialSelectParticle(target_args);

  // initialize select_mobile_
  argtype mobile_args;
  mobile_args.insert({"load_coordinates", "true"});
  mobile_args.insert({"particle_type",
    args_.key("particle_type").dflt("-1").str()});
  mobile_args.insert({"site", str(site_index_)});
  select_mobile_ = TrialSelectParticle(mobile_args);

  ASSERT(!args_.key("group_index").used(), "group not implemented with AVB");
}

class MapSelectParticleAVB {
 public:
  MapSelectParticleAVB() {
    auto obj = MakeSelectParticleAVB(MakeNeighborCriteria(),
      {{"grand_canonical", "true"}});
    obj->deserialize_map()["SelectParticleAVB"] = obj;
  }
};

static MapSelectParticleAVB mapper_ = MapSelectParticleAVB();

void SelectParticleAVB::precompute(System * system) {
  TrialSelect::precompute(system);
  select_target_.precompute(system);
  select_mobile_.precompute(system);
  anchor_.clear();
  anchor_.add_site(0, select_target_.site());
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
      config, &second_target_, random);
    if (num2 <= 0) return false;
    anchor_.set_particle(0, second_target_.particle_index(0));
  } else {
    anchor_.set_particle(0, target_.particle_index(0));
  }
  ASSERT(target_.num_sites() == 1, "Error");
  map_(*system, *neighbor_criteria_).neighbors(
    *neighbor_criteria_,
    config,
    target_.particle_index(0),
    target_.site_index(0, 0),
    site_index_,
    random,
    &neighbors_);
  const int num_neighbors = static_cast<int>(neighbors_.num_sites());
  DEBUG("neighbors: " << neighbors_.str());

  // Initialize mobile
  int num_out = 0.;
  if (!is_ghost() && inside_) {
    DEBUG("inside. num_neighbors: " << num_neighbors);
    if (num_neighbors <= 0) return false;
    mobile_.set_particle(0,
      random->const_element(neighbors_.particle_indices()));
    if (!grand_canonical_) {
      DEBUG("loading positions");
      mobile_.load_positions(config.particles());
    }
  } else if (!inside_) {
    // only select that is outside: AVB2 out->in
    DEBUG("outside");
    neighbors_.add(target_); // add target to neighbors to exclude.
    num_out = select_mobile_.random_particle(config,
      &neighbors_,
      &mobile_,
      random);
    DEBUG("num_out: " << num_out);
    if (num_out == 0) return false;
    mobile_.load_positions(config.particles());
  }

  // precompute volume terms
  const double volume = config.domain().volume();
  const double volume_av = neighbor_criteria_->volume(config.dimension());
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
      system->get_configuration(), &empty_, &mobile_);
    set_probability(volume_av/static_cast<double>(num_neighbors + 1));
    if (target_mobile_same_type_) {
      // ghost will be added during perturb. Not yet added, so n/(n+1) factor.
      set_probability(probability()
        *static_cast<double>(num)/static_cast<double>(num + 1));
    }

  // GCE remove
  } else if (!is_ghost() && grand_canonical_ && inside_ && !is_second_target_) {
    set_probability(static_cast<double>(num_neighbors)/volume_av);
    if (target_mobile_same_type_) {
      set_probability(probability()
        *static_cast<double>(num)/static_cast<double>(num - 1));
    }

  // AVB2 in->out
  } else if (!is_ghost() && !grand_canonical_ && inside_ && !is_second_target_) {
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
    set_probability(static_cast<double>(num_neighbors)
      /static_cast<double>(num_out + 1)
      *volume_out/volume_av);

  // AVB2 out->in
  } else if (!is_ghost() && !grand_canonical_ && !inside_ && !is_second_target_) {
    // ComputeAVB2 will add p_bias term
    set_probability(static_cast<double>(num_out)
      /static_cast<double>(num_neighbors + 1)
      *volume_av/volume_out);

  // AVB4 in->in
  } else if (!is_ghost() && !grand_canonical_ && inside_ && is_second_target_) {
    // obtain the number of neighbors in the second target
    map_(*system, *neighbor_criteria_).neighbors(
      *neighbor_criteria_,
      config,
      second_target_.particle_index(0),
      second_target_.site_index(0, 0),
      site_index_,
      random,
      &neighbors_);
    int mobile_not_in_both_targets_ = 1;
    if (find_in_list(mobile_.particle_index(0),
                     neighbors_.particle_indices())) {
      mobile_not_in_both_targets_ = 0;
    }
    set_probability(static_cast<double>(num_neighbors)
      /static_cast<double>(neighbors_.num_particles() +
                           mobile_not_in_both_targets_));

  // invalid
  } else {
    FATAL("unrecognized combination of the is_ghost: " << is_ghost()
      << " grand_canonical: " << grand_canonical_ << " inside: " << inside_
      << " second_target: " << is_second_target_);
  }

  ASSERT(mobile_.particle_index(0) != target_.particle_index(0),
    "mobile particle: " << mobile_.particle_index(0) <<
    " should not be equal to target: " << target_.particle_index(0));
  DEBUG("probability: " << probability());
  DEBUG("mobile: " << mobile_.str());
//  printable_["cluster_size"].accumulate(mobile_.num_particles());
  remove_unphysical_sites(config);
  mobile_original_ = mobile_;
  return true;
}

std::shared_ptr<TrialSelect> SelectParticleAVB::create(std::istream& istr) const {
  return std::make_shared<SelectParticleAVB>(istr);
}

SelectParticleAVB::SelectParticleAVB(std::istream& istr)
  : TrialSelect(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(239 == version, "mismatch version: " << version);
  // HWH for unknown reasons, this function template does not work
  // feasst_deserialize(neighbor_criteria_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      neighbor_criteria_ = std::make_shared<NeighborCriteria>(istr);
    }
  }
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
  feasst_serialize(neighbor_criteria_, ostr);
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
