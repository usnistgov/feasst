#include "utils/include/serialize.h"
#include "cluster/include/select_cluster.h"

namespace feasst {

SelectCluster::SelectCluster(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args)
  : TrialSelect(args) {
  class_name_ = "SelectCluster";
  neighbor_criteria_ = neighbor_criteria;
  select_particle_ = std::make_shared<TrialSelectParticle>(args);
  // HWH optimize by not loading coordinates?
  printable_["cluster_size"] = Accumulator();
}

class MapSelectCluster {
 public:
  MapSelectCluster() {
    auto obj = MakeSelectCluster(MakeNeighborCriteria());
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
  map_(system, neighbor_criteria_.get())->select_cluster(
    neighbor_criteria_.get(), config, first_particle,
    select, frame_of_reference);
}

bool SelectCluster::are_constraints_satisfied(const System& system) const {
  return !map_(system, neighbor_criteria_.get())->is_cluster_changed(
    neighbor_criteria_.get(), mobile_, system.configuration());
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
  set_probability(1./static_cast<double>(num));
  select_cluster(first_particle, *system);
  printable_["cluster_size"].accumulate(mobile_.num_particles());
  remove_unphysical_sites(config);
  mobile_original_ = mobile_;
  return true;
}

std::shared_ptr<TrialSelect> SelectCluster::create(std::istream& istr) const {
  return std::make_shared<SelectCluster>(istr);
}

SelectCluster::SelectCluster(std::istream& istr)
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
  feasst_serialize(neighbor_criteria_, ostr);
  feasst_serialize_fstdr(select_particle_, ostr);
}

void SelectCluster::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_select_cluster_(ostr);
}

}  // namespace feasst
