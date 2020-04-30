
#ifndef FEASST_CLUSTER_ENERGY_MAP_ALL_CRITERIA_H_
#define FEASST_CLUSTER_ENERGY_MAP_ALL_CRITERIA_H_

#include <vector>
#include "utils/include/arguments.h"
#include "system/include/neighbor_criteria.h"
#include "cluster/include/energy_map_all.h"

namespace feasst {
/**
  Same as EnergyMapAll, except subject update to NeighborCriteria.
 */
class EnergyMapAllCriteria : public EnergyMapAll {
 public:
  EnergyMapAllCriteria(std::shared_ptr<NeighborCriteria> neighbor_criteria,
      const argtype& args = argtype()) : EnergyMapAll(args) {
    neighbor_criteria_ = neighbor_criteria;
  }

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
  const NeighborCriteria * neighbor_criteria() const override {
    return neighbor_criteria_.get(); }

  // serialization
  std::string class_name() const override { return class_name_; }
  std::shared_ptr<EnergyMap> create(std::istream& istr) const override {
    return std::make_shared<EnergyMapAllCriteria>(istr); }
  void serialize(std::ostream& ostr) const override;
  EnergyMapAllCriteria(std::istream& istr);
  virtual ~EnergyMapAllCriteria() {}

 private:
  const std::string class_name_ = "EnergyMapAllCriteria";
  std::shared_ptr<NeighborCriteria> neighbor_criteria_;
};

inline std::shared_ptr<EnergyMapAllCriteria> MakeEnergyMapAllCriteria(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args = argtype()) {
  return std::make_shared<EnergyMapAllCriteria>(neighbor_criteria, args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_ENERGY_MAP_ALL_CRITERIA_H_
