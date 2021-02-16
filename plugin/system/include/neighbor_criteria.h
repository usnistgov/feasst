
#ifndef FEASST_SYSTEM_NEIGHBOR_CRITERIA_H_
#define FEASST_SYSTEM_NEIGHBOR_CRITERIA_H_

#include <memory>
#include "utils/include/arguments.h"
#include "math/include/position.h"

namespace feasst {

class Domain;

// HWH consider applying Shapes from confinement into NeighborCriteria
// For now, it assumes spheres/circles
/**
  Criteria for determining neighbors.
 */
class NeighborCriteria {
 public:
  /**
    args:
    - reference_potential: index of reference potentials (default: -1).
      If -1, use full potentials.
    - potential_index: index of potential for pair interaction (default: 0).
    - energy_maximum: maximum energy to be in cluster (default: -NEAR_ZERO).
    - minimum_distance: minimum separation distance (default: 0).
    - maximum_distance: maximum separation distance (default: NEAR_INFINITY).
    - site_type0: consider only interactions between a specific site type.
      If -1, consider all sites (default: -1).
    - site_type1: consider only interactions between a specific site type.
      If -1, consider all sites (default: -1).
   */
  explicit NeighborCriteria(argtype args = argtype());

  int reference_potential() const { return reference_potential_; }
  int potential_index() const { return potential_index_; }
  double energy_maximum() const { return energy_maximum_; }
  double minimum_distance() const;
  double maximum_distance() const;

  /// Return true if criteria are satisfied.
  bool is_accepted(const double energy,
                   const double squared_distance,
                   const int site_type0,
                   const int site_type1) const;

  /// Return the volume based on the distance criteria.
  double volume(const int dimension) const;

  /// Return true if position satisfies criteria, taking into account PBCs.
  bool is_position_accepted(
    const Position& position,
    const Domain& domain);

  /// Serialize.
  void serialize(std::ostream& ostr) const;

  /// Construct from serialization.
  explicit NeighborCriteria(std::istream& istr);

 private:
  int reference_potential_;
  int potential_index_;
  double energy_maximum_;
  double minimum_distance_sq_;
  double maximum_distance_sq_;
  int site_type0_;
  int site_type1_;

  // temporary
  Position rel_, pbc_, origin_;
};

inline std::shared_ptr<NeighborCriteria> MakeNeighborCriteria(
    argtype args = argtype()) {
  return std::make_shared<NeighborCriteria>(args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_NEIGHBOR_CRITERIA_H_
