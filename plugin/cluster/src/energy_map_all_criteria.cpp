#include "cluster/include/energy_map_all_criteria.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapEnergyMapAllCriteria {
 public:
  MapEnergyMapAllCriteria() {
    auto neighbor_criteria = MakeNeighborCriteria();
    EnergyMapAllCriteria(neighbor_criteria).deserialize_map()["EnergyMapAllCriteria"] =
      MakeEnergyMapAllCriteria(neighbor_criteria);
  }
};

static MapEnergyMapAllCriteria mapper_ = MapEnergyMapAllCriteria();

EnergyMapAllCriteria::EnergyMapAllCriteria(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args) : EnergyMapAll(args) {
  class_name_ = "EnergyMapAllCriteria";
  neighbor_criteria_ = neighbor_criteria;
}

EnergyMapAllCriteria::EnergyMapAllCriteria(std::istream& istr)
  : EnergyMapAll(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9284, "mismatch:" << version);
  // feasst_deserialize(neighbor_criteria_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      neighbor_criteria_ = std::make_shared<NeighborCriteria>(istr);
    }
  }
}

void EnergyMapAllCriteria::serialize(std::ostream& ostr) const {
  serialize_energy_map_all_(ostr);
  feasst_serialize_version(9284, ostr);
  feasst_serialize(neighbor_criteria_, ostr);
}

double EnergyMapAllCriteria::update(
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
    return EnergyMapAll::update(energy, part1_index, site1_index, site1_type,
      part2_index, site2_index, site2_type, squared_distance, pbc);
  }
  return energy;
}

}  // namespace feasst
