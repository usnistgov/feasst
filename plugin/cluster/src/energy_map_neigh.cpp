#include "cluster/include/energy_map_neigh.h"
#include "utils/include/utils.h"  // find_in_list
#include "utils/include/serialize.h"

namespace feasst {

class MapEnergyMapNeigh {
 public:
  MapEnergyMapNeigh() {
    EnergyMapNeigh().deserialize_map()["EnergyMapNeigh"] =
      MakeEnergyMapNeigh();
  }
};

static MapEnergyMapNeigh mapper_ = MapEnergyMapNeigh();

EnergyMapNeigh::EnergyMapNeigh(std::istream& istr) : EnergyMap(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 210, "mismatch:" << version);
  feasst_deserialize(&map_, istr);
}

void EnergyMapNeigh::serialize_energy_map_neigh_(std::ostream& ostr) const {
  feasst_serialize_version(210, ostr);
  feasst_serialize(map_, ostr);
}

void EnergyMapNeigh::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_energy_map_neigh_(ostr);
}

std::vector<double> * EnergyMapNeigh::smap_(const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index) {
  int index;
  const bool found = find_in_list(part2_index, neighbor_[part1_index], &index);
  ASSERT(found, "not found");
  return &map_[part1_index][index][site1_index][site2_index];
}

std::vector<double> * EnergyMapNeigh::smap_new_(const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index) {
  int index;
  const bool found = find_in_list(part2_index, neighbor_new_[part1_index], &index);
  ASSERT(found, "not found");
  return &map_new_[part1_index][index][site1_index][site2_index];
}

double EnergyMapNeigh::update(
    const double energy,
    const int part1_index,
    const int site1_index,
    const int site1_type,
    const int part2_index,
    const int site2_index,
    const int site2_type,
    const double squared_distance,
    const Position * pbc) {
  FATAL("not implemented");
}

//void EnergyMapNeigh::remove_particles(const Select& select) {
//  DEBUG("sel: " << select.str());
//  const int pmax = select.particle_indices().back();
//  const int smax = select.site_indices().back().back();
//  resize_(pmax, smax, pmax, smax);
//  DEBUG("map size: " << map_.size());
//  DEBUG("map size: " << map_[0].size());
//  DEBUG("map size: " << map_[0][0].size());
//  DEBUG("map size: " << map_[0][0][0].size());
//  DEBUG("map: " << feasst_str(map_));
//  for (int sel1 = 0; sel1 < select.num_particles(); ++sel1) {
//    const int p1 = select.particle_index(sel1);
//    // For each other particle, remove p1 as neighbor to others
//    // But for neighbor_[p1], simply empty it and leave blank
//    for (int p2 : neighbor_[p1]) {
//      // find index where neighbor_[p2][index] == p1
//      // remove neighbor_[p2][index]
//    }
//    // clear neighbor_[p1]
//  }
//}

void EnergyMapNeigh::revert(const Select& select) {
  for (int sel_index = 0; sel_index < select.num_particles(); ++sel_index) {
    const int p1 = select.particle_index(sel_index);
    for (const int s1 : select.site_indices(sel_index)) {
      for (int p2 = 0; p2 < static_cast<int>(map_[p1].size()); ++p2) {
        for (int s2 = 0; s2 < static_cast<int>(map_[p1][p2][s1].size()); ++s2) {
          map_new_[p1][p2][s1][s2] = map_[p1][p2][s1][s2];
          map_new_[p2][p1][s2][s1] = map_[p2][p1][s2][s1];
        }
      }
    }
  }
}

void EnergyMapNeigh::select_cluster(
    const NeighborCriteria * neighbor_criteria,
    const Configuration& config,
    const int particle_node,
    Select * cluster,
    const Position& frame_of_reference) const {
  FATAL("not impl");
}

bool EnergyMapNeigh::is_cluster_(
    const NeighborCriteria * neighbor_criteria,
    const std::vector<std::vector<std::vector<double> > >& smap,
    const Configuration& config,
    Position * frame) const {
  FATAL("not impl");
  return false;
}

void EnergyMapNeigh::resize_(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index) { FATAL("not impl"); }

}  // namespace feasst
