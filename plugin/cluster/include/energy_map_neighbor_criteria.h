
#ifndef FEASST_CLUSTER_ENERGY_MAP_NEIGHBOR_CRITERIA_H_
#define FEASST_CLUSTER_ENERGY_MAP_NEIGHBOR_CRITERIA_H_

#include <vector>
#include "utils/include/arguments.h"
#include "system/include/neighbor_criteria.h"
#include "cluster/include/energy_map_neighbor.h"

namespace feasst {
/**
  Same as EnergyMapNeighbor, except subject update to NeighborCriteria.
 */
class EnergyMapNeighborCriteria : public EnergyMapNeighbor {
 public:
  EnergyMapNeighborCriteria(std::shared_ptr<NeighborCriteria> neighbor_criteria,
      const argtype& args = argtype());

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
  bool is_queryable() const override { return false; }
  const NeighborCriteria& neighbor_criteria() const override {
    return *neighbor_criteria_; }

  // serialization
  std::string class_name() const override { return class_name_; }
  std::shared_ptr<EnergyMap> create(std::istream& istr) const override {
    return std::make_shared<EnergyMapNeighborCriteria>(istr); }
  void serialize(std::ostream& ostr) const override;
  EnergyMapNeighborCriteria(std::istream& istr);
  virtual ~EnergyMapNeighborCriteria() {}

 private:
  std::shared_ptr<NeighborCriteria> neighbor_criteria_;
};

inline std::shared_ptr<EnergyMapNeighborCriteria> MakeEnergyMapNeighborCriteria(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args = argtype()) {
  return std::make_shared<EnergyMapNeighborCriteria>(neighbor_criteria, args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_ENERGY_MAP_NEIGHBOR_CRITERIA_H_
