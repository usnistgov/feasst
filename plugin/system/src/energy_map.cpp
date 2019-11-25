#include "system/include/energy_map.h"

namespace feasst {

class MapEnergyMap {
 public:
  MapEnergyMap() {
    EnergyMap().deserialize_map()["EnergyMap"] =
      std::make_shared<EnergyMap>();
  }
};

static MapEnergyMap mapper_energy_map_ = MapEnergyMap();

std::map<std::string, std::shared_ptr<EnergyMap> >& EnergyMap::deserialize_map() {
  static std::map<std::string, std::shared_ptr<EnergyMap> >* ans =
     new std::map<std::string, std::shared_ptr<EnergyMap> >();
  return *ans;
}

//
//
//

class MapEnergyMapAll {
 public:
  MapEnergyMapAll() {
    EnergyMapAll().deserialize_map()["EnergyMapAll"] =
      std::make_shared<EnergyMapAll>();
  }
};

static MapEnergyMapAll mapper_ = MapEnergyMapAll();

EnergyMapAll::EnergyMapAll(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 210, "mismatch:" << version);
  feasst_deserialize(&map_, istr);
  feasst_deserialize(&site_max_, istr);
}

void EnergyMapAll::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(210, ostr);
  feasst_serialize(map_, ostr);
  feasst_serialize(site_max_, ostr);
}

void EnergyMapAll::precompute(Configuration * config) {
  site_max_ = config->max_sites_in_any_particle();
  DEBUG("site_max_ " << site_max_);
}

double EnergyMapAll::update(
    const double energy,
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index) {
  // resize and initialize
  const int part_max_old = static_cast<int>(map_.size());
  if (part1_index >= part_max_old || part2_index >= part_max_old) {
    const int part_max = std::max(part1_index, part2_index);
    map_.resize(part_max + 1);
    for (int i = 0; i < part_max + 1; ++i) {
      map_[i].resize(part_max + 1);
      for (int j = 0; j < part_max + 1; ++j) {
        ASSERT(site_max_ != 0, "wasn't precomputed");
        map_[i][j].resize(site_max_ + 1);
        for (int k = 0; k < site_max_ + 1; ++k) {
          map_[i][j][k].resize(site_max_ + 1);
          if (i > part_max_old || j > part_max_old) {
            for (int l = 0; l < site_max_ + 1; ++l) {
              map_[i][j][k][l] = 0.;
            }
          }
        }
      }
    }
  }
//  INFO(part1_index << " " << part2_index << " " << site1_index << " " << site2_index);
//  INFO(map_.size() << " " << map_[0].size() << " " << map_[0][0].size() << " " << map_[0][0][0].size());
  map_[part1_index][part2_index][site1_index][site2_index] = energy;
  map_[part2_index][part1_index][site2_index][site1_index] = energy;
  return energy;
}

//void EnergyMapAll::revert() {
//  INFO("hi");
//}

}  // namespace feasst
