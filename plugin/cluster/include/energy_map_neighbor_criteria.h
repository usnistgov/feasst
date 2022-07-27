
#ifndef FEASST_CLUSTER_ENERGY_MAP_NEIGHBOR_CRITERIA_H_
#define FEASST_CLUSTER_ENERGY_MAP_NEIGHBOR_CRITERIA_H_

#include <vector>
#include "utils/include/arguments.h"
#include "configuration/include/neighbor_criteria.h"
#include "cluster/include/energy_map_neighbor.h"

namespace feasst {
/**
  Same as EnergyMapNeighbor, except subject update to NeighborCriteria.
 */
class EnergyMapNeighborCriteria : public EnergyMapNeighbor {
 public:
  /**
    args:
    - neighbor_index: NeighborCriteria index contained in Configuration (default: 0).
   */
  explicit EnergyMapNeighborCriteria(argtype args = argtype());
  explicit EnergyMapNeighborCriteria(argtype * args);

  double update(
      const double energy,
      const int part1_index,
      const int site1_index,
      const int site1_type,
      const int part2_index,
      const int site2_index,
      const int site2_type,
      const double squared_distance,
      const Position * pbc,
      const Configuration& config) override;
  bool is_queryable() const override { return false; }

  // serialization
  std::string class_name() const override { return class_name_; }
  std::shared_ptr<EnergyMap> create(std::istream& istr) const override {
    return std::make_shared<EnergyMapNeighborCriteria>(istr); }
  std::shared_ptr<EnergyMap> create(argtype * args) const override {
    return std::make_shared<EnergyMapNeighborCriteria>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit EnergyMapNeighborCriteria(std::istream& istr);
  virtual ~EnergyMapNeighborCriteria() {}

 private:
  int neighbor_index_;
};

inline std::shared_ptr<EnergyMapNeighborCriteria> MakeEnergyMapNeighborCriteria(
  argtype args = argtype()) {
  return std::make_shared<EnergyMapNeighborCriteria>(args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_ENERGY_MAP_NEIGHBOR_CRITERIA_H_
