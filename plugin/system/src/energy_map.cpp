#include "system/include/energy_map.h"

namespace feasst {

std::map<std::string, std::shared_ptr<EnergyMap> >& EnergyMap::deserialize_map() {
  static std::map<std::string, std::shared_ptr<EnergyMap> >* ans =
     new std::map<std::string, std::shared_ptr<EnergyMap> >();
  return *ans;
}

}  // namespace feasst
