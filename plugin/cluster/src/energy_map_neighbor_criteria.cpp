#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "configuration/include/configuration.h"
#include "cluster/include/energy_map_neighbor_criteria.h"

namespace feasst {

FEASST_MAPPER(EnergyMapNeighborCriteria,);

EnergyMapNeighborCriteria::EnergyMapNeighborCriteria(argtype * args)
  : EnergyMapNeighbor(args) {
  class_name_ = "EnergyMapNeighborCriteria";
  neighbor_index_ = integer("neighbor_index", args, 0);
}
EnergyMapNeighborCriteria::EnergyMapNeighborCriteria(argtype args)
  : EnergyMapNeighborCriteria(&args) {
  feasst_check_all_used(args);
}

EnergyMapNeighborCriteria::EnergyMapNeighborCriteria(std::istream& istr)
  : EnergyMapNeighbor(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3466, "mismatch:" << version);
  feasst_deserialize(&neighbor_index_, istr);
}

void EnergyMapNeighborCriteria::serialize(std::ostream& ostr) const {
  serialize_energy_map_neighbor_(ostr);
  feasst_serialize_version(3466, ostr);
  feasst_serialize(neighbor_index_, ostr);
}

double EnergyMapNeighborCriteria::update(
    const double energy,
    const int part1_index,
    const int site1_index,
    const int site1_type,
    const int part2_index,
    const int site2_index,
    const int site2_type,
    const double squared_distance,
    const Position * pbc,
    const Configuration& config) {
  TRACE("updating neighbor_index:" << neighbor_index_);
  if (config.neighbor_criteria(neighbor_index_).is_accepted(energy, squared_distance,
                                      site1_type, site2_type)) {
    return EnergyMapNeighbor::update(energy, part1_index, site1_index, site1_type,
      part2_index, site2_index, site2_type, squared_distance, pbc, config);
  }
  finalizable_ = true;
  return energy;
}

}  // namespace feasst
