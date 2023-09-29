#include "utils/include/serialize.h"
#include "cluster/include/calculate_cluster.h"

namespace feasst {

class MapCalculateCluster {
 public:
  MapCalculateCluster() {
    auto obj = MakeCalculateCluster();
    obj->deserialize_map()["CalculateCluster"] = obj;
  }
};

static MapCalculateCluster mapper_ = MapCalculateCluster();

CalculateCluster::CalculateCluster(argtype * args) : Modify(args) {}
CalculateCluster::CalculateCluster(argtype args) : CalculateCluster(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void CalculateCluster::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  printer(header(*criteria, *system, *trial_factory),
          file_name(*criteria));
}

std::string CalculateCluster::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  ss << accumulator_.status_header() << std::endl;
  return ss.str();
}

void CalculateCluster::update(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  const Configuration& config = system->configuration();
  if (config.num_particles() > 0) {
    INFO(config.selection_of_all().particle_index(0));
  }
  const int neigh = 0;
  const int conf = 0;
  Select sel_all;
  std::vector<Select> cluster;
  system->energy();
  for (int first_particle : config.selection_of_all().particle_indices()) {
    Select select;
    const Particle& part = config.select_particle(first_particle);
    if (!find_in_list(first_particle, sel_all.particle_indices())) {
    //if (!sel_all.is_overlap(select)) {
      INFO("no overlap?");
      Position frame_of_reference(config.dimension());
      select.add_particle(part, first_particle);
      select.load_positions_of_last(part, frame_of_reference);
      system->potential(0).visit_model().inner().energy_map().select_cluster(
        system->neighbor_criteria(neigh, conf),
        system->configuration(conf),
        first_particle,
        &select,
        frame_of_reference
      );
      INFO(select.str());
      accumulator_.accumulate(select.num_particles());
      sel_all.add(select);
      cluster.push_back(select);
    }
  }
  INFO(sel_all.num_particles());
  ASSERT(sel_all.num_particles() == config.num_particles(),
    sel_all.num_particles() << " != " << config.num_particles());
}

std::string CalculateCluster::write(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(*criteria, *system, *trial_factory);
  }
  ss << accumulator_.status() << std::endl;
  DEBUG(ss.str());
  return ss.str();
}

void CalculateCluster::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(8349, ostr);
}

CalculateCluster::CalculateCluster(std::istream& istr) : Modify(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8349, "mismatch version:" << version);
}

CalculateCluster::CalculateCluster(const Modify& energy) {
  std::stringstream ss;
  energy.serialize(ss);
  *this = CalculateCluster(ss);
}

}  // namespace feasst
