#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "configuration/include/configuration.h"
#include "cluster/include/energy_map_all_criteria.h"

namespace feasst {

class MapEnergyMapAllCriteria {
 public:
  MapEnergyMapAllCriteria() {
    auto obj = MakeEnergyMapAllCriteria();
    obj->deserialize_map()["EnergyMapAllCriteria"] = obj;
  }
};

static MapEnergyMapAllCriteria mapper_ = MapEnergyMapAllCriteria();

EnergyMapAllCriteria::EnergyMapAllCriteria(argtype * args) : EnergyMapAll(args) {
  class_name_ = "EnergyMapAllCriteria";
  neighbor_index_ = integer("neighbor_index", args, 0);
}
EnergyMapAllCriteria::EnergyMapAllCriteria(argtype args) : EnergyMapAllCriteria(&args) {
  feasst_check_all_used(args);
}

EnergyMapAllCriteria::EnergyMapAllCriteria(std::istream& istr)
  : EnergyMapAll(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9284, "mismatch:" << version);
  feasst_deserialize(&neighbor_index_, istr);
}

void EnergyMapAllCriteria::serialize(std::ostream& ostr) const {
  serialize_energy_map_all_(ostr);
  feasst_serialize_version(9284, ostr);
  feasst_serialize(neighbor_index_, ostr);
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
    const Position * pbc,
    const Configuration& config) {
  if (config.neighbor_criteria(neighbor_index_).is_accepted(energy,
    squared_distance, site1_type, site2_type)) {
    return EnergyMapAll::update(energy, part1_index, site1_index, site1_type,
      part2_index, site2_index, site2_type, squared_distance, pbc, config);
  }
  return energy;
}

}  // namespace feasst
