#include "utils/include/serialize.h"
#include "cluster/include/select_cluster.h"

namespace feasst {

SelectCluster::SelectCluster(argtype args) : SelectCluster(&args) {
  FEASST_CHECK_ALL_USED(args);
}
SelectCluster::SelectCluster(argtype * args) : TrialSelect(args) {
  class_name_ = "SelectCluster";
  neighbor_ = integer("neighbor_index", args, 0);
  if (is_particle_type_set()) {
    args->insert({"particle_type", str(particle_type())});
  }
  select_particle_ = std::make_shared<TrialSelectParticle>(args);
  printable_["cluster_size"] = Accumulator();
}

class MapSelectCluster {
 public:
  MapSelectCluster() {
    auto obj = MakeSelectCluster();
    obj->deserialize_map()["SelectCluster"] = obj;
  }
};

static MapSelectCluster mapper_ = MapSelectCluster();

void SelectCluster::precompute(System * system) {
  TrialSelect::precompute(system);
  select_particle_->precompute(system);
}

void SelectCluster::select_cluster(const int first_particle,
    const System& system,
    Select * select) {
  select->clear();
  const Configuration& config = system.configuration();
  const Particle& part = config.select_particle(first_particle);
  Position frame_of_reference(config.dimension());
  select->add_particle(part, first_particle);
  select->load_positions_of_last(part, frame_of_reference);
  DEBUG("first node " << select->str());
  DEBUG("first node pos " << select->site_positions()[0][0].str());
  // HWH update with configuration_index_
  map_(system, neighbor_).select_cluster(system.neighbor_criteria(neighbor_, 0),
    config, first_particle, select, frame_of_reference);
}

bool SelectCluster::are_constraints_satisfied(const int old,
    const System& system) const {
  if (old == 0) {
    // HWH update with configuration_index_
    const bool constraint = !map_(system, neighbor_).is_cluster_changed(
      system.neighbor_criteria(neighbor_, 0), mobile_, system.configuration());
    DEBUG("constraint " << constraint);
    return constraint;
  }
  return true;
}

std::vector<Select> SelectCluster::select_clusters(
    const System& system) {
  std::vector<Select> clusters;
  Select selectable = system.configuration().group_select(group_index());
  int num_iter = 0;
  const int max_iter = selectable.num_particles();
  while (selectable.num_particles() != 0) {
    clusters.push_back(Select());
    select_cluster(selectable.particle_index(0), system, &clusters.back());
    selectable.remove(clusters.back());
    ++num_iter;
    ASSERT(num_iter <= max_iter, "maximum iterations reached");
  }
  return clusters;
}

bool SelectCluster::select(const Select& perturbed,
                           System * system,
                           Random * random) {
  const Configuration& config = system->configuration();
  Select first_node;
  // HWH use TrialSelectParticle for optional group, particle type, etc.
  const int num = select_particle_->random_particle(config, &first_node, random);
  if (num <= 0) return false;
  const int first_particle = first_node.particle_index(0);
  set_probability_(1./static_cast<double>(num));
  select_cluster(first_particle, *system);
  printable_["cluster_size"].accumulate(mobile_.num_particles());
  if (mobile_.num_particles() == 1) {
    return false;
  }
  remove_unphysical_sites(config);
  set_mobile_original(system);
  return true;
}

std::shared_ptr<TrialSelect> SelectCluster::create(std::istream& istr) const {
  return std::make_shared<SelectCluster>(istr);
}

SelectCluster::SelectCluster(std::istream& istr)
  : TrialSelect(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(239 == version, "mismatch version: " << version);
  feasst_deserialize(&neighbor_, istr);
  // HWH for unknown reasons, this function template does not work
  // feasst_deserialize_fstdr(select_particle_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      select_particle_ = std::make_shared<TrialSelectParticle>(istr);
    }
  }
}

void SelectCluster::serialize_select_cluster_(std::ostream& ostr) const {
  serialize_trial_select_(ostr);
  feasst_serialize_version(239, ostr);
  feasst_serialize(neighbor_, ostr);
  feasst_serialize_fstdr(select_particle_, ostr);
}

void SelectCluster::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_select_cluster_(ostr);
}

}  // namespace feasst
