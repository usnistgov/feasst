
#ifndef FEASST_CLUSTER_ENERGY_MAP_ALL_CRITERIA_H_
#define FEASST_CLUSTER_ENERGY_MAP_ALL_CRITERIA_H_

#include <vector>
#include "configuration/include/neighbor_criteria.h"
#include "cluster/include/energy_map_all.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Same as EnergyMapAll, except subject update to NeighborCriteria.
 */
class EnergyMapAllCriteria : public EnergyMapAll {
 public:
  //@{
  /** @name Arguments
    - neighbor_index: NeighborCriteria index contained in Configuration (default: 0).
   */
  explicit EnergyMapAllCriteria(argtype args = argtype());
  explicit EnergyMapAllCriteria(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

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
    return std::make_shared<EnergyMapAllCriteria>(istr); }
  std::shared_ptr<EnergyMap> create(argtype * args) const override {
    return std::make_shared<EnergyMapAllCriteria>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit EnergyMapAllCriteria(std::istream& istr);
  virtual ~EnergyMapAllCriteria() {}

  //@}
 private:
  int neighbor_index_;
};

inline std::shared_ptr<EnergyMapAllCriteria> MakeEnergyMapAllCriteria(
    argtype args = argtype()) {
  return std::make_shared<EnergyMapAllCriteria>(args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_ENERGY_MAP_ALL_CRITERIA_H_
