#include "cluster/include/energy_map_all_criteria.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapEnergyMapAllCriteria {
 public:
  MapEnergyMapAllCriteria() {
    auto cluster_criteria = MakeClusterCriteria();
    EnergyMapAllCriteria(cluster_criteria).deserialize_map()["EnergyMapAllCriteria"] =
      MakeEnergyMapAllCriteria(cluster_criteria);
  }
};

static MapEnergyMapAllCriteria mapper_ = MapEnergyMapAllCriteria();

EnergyMapAllCriteria::EnergyMapAllCriteria(std::istream& istr)
  : EnergyMapAll(istr) {
  // feasst_deserialize(cluster_criteria_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      cluster_criteria_ = std::make_shared<ClusterCriteria>(istr);
    }
  }
}

void EnergyMapAllCriteria::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_energy_map_all_(ostr);
  feasst_serialize(cluster_criteria_, ostr);
}

double EnergyMapAllCriteria::update(
    const double energy,
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const double squared_distance,
    const Position * pbc) {
  if (cluster_criteria_->is_accepted(energy, squared_distance)) {
    return EnergyMapAll::update(energy, part1_index, site1_index, part2_index,
                                site2_index, squared_distance, pbc);
  }
  return energy;
}

}  // namespace feasst
