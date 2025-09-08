
#ifndef FEASST_CONFIGURATION_NEIGHBOR_CRITERIA_H_
#define FEASST_CONFIGURATION_NEIGHBOR_CRITERIA_H_

#include <map>
#include <memory>
#include <string>

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class Domain;
class ParticleFactory;
class Position;

// HWH consider applying Shapes from confinement into NeighborCriteria
// For now, it assumes spheres/circles
// HWH Implement a NeighborCriteriaFactory to allow for AV definitions of
// multiple site types with different distances. OR,
// HWH Implement a csv for energy, dist and type (used in is_accepted).
// But no, each Trial needs to know volume for acceptance
// What I need to support is multiple NeighborCriteria in EnergyMap ?
/**
  Criteria for defining neighbors.
 */
class NeighborCriteria {
 public:
  //@{
  /** @name Arguments
    - ref: name of RefPotential. If empty, use Potential (default: empty)
    - potential_index: index of potential for pair interaction (default: 0).
    - energy_maximum: maximum energy to be in cluster (default: largest double precision).
    - minimum_distance: minimum separation distance (default: 0).
    - maximum_distance: maximum separation distance (default: NEAR_INFINITY).
    - site_type0: consider only interactions between a specific site type.
      If -1, consider all sites (default: -1).
      Otherwise, site_type1 must also be included.
    - site_type1: consider only interactions between a specific site type.
      If -1, consider all sites (default: -1).
    - site_type0_alt: consider interactions between another specific site type.
      If -1, ignore (default: -1).
      Otherwise, site_type1_alt must also be included.
    - site_type1_alt: consider interactions between another specific site type.
      If -1, ignore (default: -1).
   */
  explicit NeighborCriteria(argtype args = argtype());
  explicit NeighborCriteria(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Convert site type names to intergers for optimization.
  /// Otherwise, assumes name is index (and index is -1 if stoi fails).
  void name_to_index(const ParticleFactory& unique_types);

  int reference_potential() const;
  const std::string& ref() const { return ref_; }
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

  //@}
 private:
  int reference_potential_, potential_index_;
  std::string ref_;
  double energy_maximum_, minimum_distance_sq_, maximum_distance_sq_;
  int site_type0_, site_type1_, site_type0_alt_, site_type1_alt_;
  std::string site_type0_name_, site_type1_name_;
  std::string site_type0_alt_name_, site_type1_alt_name_;

  // temporary
  std::shared_ptr<Position> rel_, pbc_, origin_;
};

inline std::shared_ptr<NeighborCriteria> MakeNeighborCriteria(
    argtype args = argtype()) {
  return std::make_shared<NeighborCriteria>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_NEIGHBOR_CRITERIA_H_
