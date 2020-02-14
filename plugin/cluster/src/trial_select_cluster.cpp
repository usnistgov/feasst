#include "cluster/include/trial_select_cluster.h"

namespace feasst {

TrialSelectCluster::TrialSelectCluster(
    std::shared_ptr<ClusterCriteria> cluster_criteria,
    const argtype& args)
  : TrialSelect(args) {
  class_name_ = "TrialSelectCluster";
  cluster_criteria_ = cluster_criteria;
  select_particle_ = std::make_shared<TrialSelectParticle>(args);
  // HWH optimize by not loading coordinates?
  printable_["cluster_size"] = Accumulator();
}

class MapTrialSelectCluster {
 public:
  MapTrialSelectCluster() {
    auto obj = MakeTrialSelectCluster(MakeClusterCriteria());
    obj->deserialize_map()["TrialSelectCluster"] = obj;
  }
};

static MapTrialSelectCluster mapper_ = MapTrialSelectCluster();

const EnergyMap * TrialSelectCluster::map_(const System& system) const {
  if (cluster_criteria_->reference_potential() == -1) {
    return system.const_potentials()->potentials()[
      cluster_criteria_->potential_index()].visit_model()->inner()->energy_map();
  }
  return system.reference(cluster_criteria_->reference_potential(),
                          cluster_criteria_->potential_index()
                         ).visit_model()->inner()->energy_map();
}

void TrialSelectCluster::select_cluster(const int first_particle,
    const System& system,
    SelectList * select) {
  select->clear();
  const Configuration& config = system.configuration();
  const Particle& part = config.select_particle(first_particle);
  Position frame_of_reference(config.dimension());
  select->add_particle(part, first_particle);
  select->load_positions_of_last(part, frame_of_reference);
  DEBUG("first node " << select->str());
  DEBUG("first node pos " << select->site_positions()[0][0].str());
  map_(system)->select_cluster(cluster_criteria_.get(), config, first_particle,
                               select, frame_of_reference);
}

bool TrialSelectCluster::are_constraints_satisfied(const System& system) const {
  return !map_(system)->is_cluster_changed(cluster_criteria_.get(), mobile_);
}

std::vector<SelectList> TrialSelectCluster::select_clusters(
    const System& system) {
  std::vector<SelectList> clusters;
  SelectGroup selectable = system.configuration().group_select(group_index());
  int num_iter = 0;
  const int max_iter = selectable.num_particles();
  while (selectable.num_particles() != 0) {
    clusters.push_back(SelectList());
    select_cluster(selectable.particle_index(0), system, &clusters.back());
    selectable.remove(clusters.back());
    ++num_iter;
    ASSERT(num_iter <= max_iter, "maximum iterations reached");
  }
  return clusters;
}

bool TrialSelectCluster::select(const Select& perturbed,
                                 System* system,
                                 Random * random) {
  const Configuration& config = system->configuration();
  SelectList first_node;
  // HWH use TrialSelectParticle for optional group, particle type, etc.
  const int num = select_particle_->random_particle(config, &first_node, random);
  const int first_particle = first_node.particle_index(0);
  if (num <= 0) return false;
  set_probability(1./static_cast<double>(num));
  select_cluster(first_particle, *system);
  printable_["cluster_size"].accumulate(mobile_.num_particles());
  mobile_.remove_unphysical_sites(config);
  mobile_original_ = mobile_;
  return true;
}

std::shared_ptr<TrialSelect> TrialSelectCluster::create(std::istream& istr) const {
  return std::make_shared<TrialSelectCluster>(istr);
}

TrialSelectCluster::TrialSelectCluster(std::istream& istr)
  : TrialSelect(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(239 == version, "mismatch version: " << version);
  // HWH for unknown reasons, this function template does not work
  // feasst_deserialize(cluster_criteria_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      cluster_criteria_ = std::make_shared<ClusterCriteria>(istr);
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

void TrialSelectCluster::serialize_trial_select_cluster_(std::ostream& ostr) const {
  serialize_trial_select_(ostr);
  feasst_serialize_version(239, ostr);
  feasst_serialize(cluster_criteria_, ostr);
  feasst_serialize_fstdr(select_particle_, ostr);
}

void TrialSelectCluster::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_cluster_(ostr);
}

}  // namespace feasst
