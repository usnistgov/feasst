#include "cluster/include/energy_map_neighbor_criteria.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapEnergyMapNeighborCriteria {
 public:
  MapEnergyMapNeighborCriteria() {
    auto neighbor_criteria = MakeNeighborCriteria();
    EnergyMapNeighborCriteria(neighbor_criteria).deserialize_map()["EnergyMapNeighborCriteria"] =
      MakeEnergyMapNeighborCriteria(neighbor_criteria);
  }
};

static MapEnergyMapNeighborCriteria mapper_ = MapEnergyMapNeighborCriteria();

EnergyMapNeighborCriteria::EnergyMapNeighborCriteria(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args) : EnergyMapNeighbor(args) {
  class_name_ = "EnergyMapNeighborCriteria";
  neighbor_criteria_ = neighbor_criteria;
}

EnergyMapNeighborCriteria::EnergyMapNeighborCriteria(std::istream& istr)
  : EnergyMapNeighbor(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3466, "mismatch:" << version);
  // feasst_deserialize(neighbor_criteria_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      neighbor_criteria_ = std::make_shared<NeighborCriteria>(istr);
    }
  }
}

void EnergyMapNeighborCriteria::serialize(std::ostream& ostr) const {
  serialize_energy_map_neighbor_(ostr);
  feasst_serialize_version(3466, ostr);
  feasst_serialize(neighbor_criteria_, ostr);
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
    const Position * pbc) {
  if (neighbor_criteria_->is_accepted(energy, squared_distance,
                                      site1_type, site2_type)) {
    return EnergyMapNeighbor::update(energy, part1_index, site1_index, site1_type,
      part2_index, site2_index, site2_type, squared_distance, pbc);
  }
  return energy;
}

}  // namespace feasst
