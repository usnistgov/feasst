#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "utils/include/utils.h"
#include "configuration/include/select.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "system/include/visit_model.h"
#include "system/include/potential.h"
#include "system/include/visit_model_inner.h"
#include "system/include/energy_map.h"
#include "monte_carlo/include/monte_carlo.h"
#include "cluster/include/calculate_cluster.h"

namespace feasst {

FEASST_MAPPER(CalculateCluster,);

CalculateCluster::CalculateCluster(argtype * args) : Modify(args) {}
CalculateCluster::CalculateCluster(argtype args) : CalculateCluster(&args) {
  feasst_check_all_used(args);
}

void CalculateCluster::initialize(MonteCarlo * mc) {
  Modify::initialize(mc);
  printer(header(*mc), output_file(mc->criteria()));
}

std::string CalculateCluster::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  ss << accumulator().status_header() << std::endl;
  return ss.str();
}

void CalculateCluster::update(MonteCarlo * mc) {
  System * system = mc->get_system();
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
      get_accumulator()->accumulate(select.num_particles());
      sel_all.add(select);
      cluster.push_back(select);
    }
  }
  INFO(sel_all.num_particles());
  ASSERT(sel_all.num_particles() == config.num_particles(),
    sel_all.num_particles() << " != " << config.num_particles());
}

std::string CalculateCluster::write(MonteCarlo * mc) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(*mc);
  }
  ss << accumulator().status() << std::endl;
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
