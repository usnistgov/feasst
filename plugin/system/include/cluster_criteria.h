
#ifndef FEASST_SYSTEM_CLUSTER_CRITERIA_H_
#define FEASST_SYSTEM_CLUSTER_CRITERIA_H_

#include "utils/include/arguments.h"
#include "system/include/model.h"
#include "configuration/include/configuration.h"
#include "system/include/select_list.h"

namespace feasst {

/**
  Criteria for determining clusters.
 */
class ClusterCriteria {
 public:
  /**
    args:
    - reference_potential: index of reference potentials (default: -1).
      If -1, use full potentials.
    - potential_index: index of potential for pair interaction (default: 0).
    - energy_maximum: minimum energy to be in cluster (default: 0).
    - minimum_distance: minimum separation distance (default: 0).
    - maximum_distance: maximum separation distance (default: NEAR_INFINITY).
   */
  ClusterCriteria(const argtype& args = argtype());

  int reference_potential() const { return reference_potential_; }
  int potential_index() const { return potential_index_; }
  double energy_maximum() const { return energy_maximum_; }
  double minimum_distance() const { return std::sqrt(minimum_distance_sq_); }
  double maximum_distance() const { return std::sqrt(maximum_distance_sq_); }

  /// Return true if criteria are satisfied.
  bool is_accepted(
    /// data is assumed to be in the format of EnergyMap
    /// (e.g., energy, squared_distance, pbc shifts).
    const std::vector<double>& data) const;

  /// Serialize.
  void serialize(std::ostream& ostr) const;

  /// Construct from serialization.
  ClusterCriteria(std::istream& istr);

 private:
  int reference_potential_;
  int potential_index_;
  double energy_maximum_;
  double minimum_distance_sq_;
  double maximum_distance_sq_;
};

inline std::shared_ptr<ClusterCriteria> MakeClusterCriteria(
    const argtype& args = argtype()) {
  return std::make_shared<ClusterCriteria>(args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_CLUSTER_CRITERIA_H_
