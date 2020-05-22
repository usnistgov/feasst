
#ifndef FEASST_CLUSTER_ENERGY_MAP_NEIGHBOR_H_
#define FEASST_CLUSTER_ENERGY_MAP_NEIGHBOR_H_

#include <vector>
#include "utils/include/arguments.h"
#include "system/include/energy_map.h"
#include "system/include/neighbor_criteria.h"

namespace feasst {
/**
  Map energies between neighboring particles.
  This data structure is intended to scale with particles
  better than EnergyMapAll.
  Although for small system sizes or large cutoffs, EnergyMapAll may be faster
  because it does not require sorting.

  Clear is used to generate a new map with only those particles in part1
  that represent the computed selection.

  Update populates the new map, which should be appropriately sized for part1.

  Finalize uses the new map to update the map.
 */
class EnergyMapNeighbor : public EnergyMap {
 public:
  EnergyMapNeighbor(const argtype& args = argtype());
  double update(
      const double energy,
      const int part1_index,
      const int site1_index,
      const int site1_type,
      const int part2_index,
      const int site2_index,
      const int site2_type,
      const double squared_distance,
      const Position * pbc) override;
  void revert(const Select& select) override;
  void select_cluster(const NeighborCriteria& neighbor_criteria,
                      const Configuration& config,
                      const int particle_node,
                      Select * cluster,
                      const Position& frame_of_reference) const override;
  void clear(
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index) override {}

  // serialization
  std::string class_name() const override { return class_name_; }
  std::shared_ptr<EnergyMap> create(std::istream& istr) const override {
    return std::make_shared<EnergyMapNeighbor>(istr); }
  void serialize(std::ostream& ostr) const override;
  EnergyMapNeighbor(std::istream& istr);
  virtual ~EnergyMapNeighbor() {}

 protected:
  void serialize_energy_map_neighbor_(std::ostream& ostr) const;
  void resize_(const int part1, const int site1, const int part2, const int site2) override;
  std::vector<double> * smap_(const int part1_index,
                              const int site1_index,
                              const int part2_index,
                              const int site2_index) override;
  std::vector<double> * smap_new_(const int part1_index,
                                  const int site1_index,
                                  const int part2_index,
                                  const int site2_index) override;
  const std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > >& map() const override { return map_; }

 private:
  /// The first index is the particle index, which mirrors config
  /// The second index for the list of neighbors
  /// The data are the particle index from config
  /// The neighbors should be sorted to allow for quick comparison with
  /// other lists of neighbors.
  std::vector<std::vector<int> > neighbor_, neighbor_new_;

  /// As opposed to EnergyMapAll, the second index corresponds with neighbor_ above,
  /// instead of listing all particles in config
  std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > map_, map_new_;

  int part_max_() { return static_cast<int>(map_.size()); }
  bool is_cluster_(const NeighborCriteria& neighbor_criteria,
                   const std::vector<std::vector<std::vector<double> > >& smap,
                   const Configuration& config,
                   Position * frame) const;
};

inline std::shared_ptr<EnergyMapNeighbor> MakeEnergyMapNeighbor(
    const argtype& args = argtype()) {
  return std::make_shared<EnergyMapNeighbor>(args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_ENERGY_MAP_NEIGHBOR_H_
