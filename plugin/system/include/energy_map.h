
#ifndef FEASST_SYSTEM_ENERGY_MAP_H_
#define FEASST_SYSTEM_ENERGY_MAP_H_

#include <string>
#include <vector>
#include <memory>
#include <map>
#include "system/include/synchronize_data.h"

namespace feasst {

class Random;
class Position;
class Configuration;
class Select;
class NeighborCriteria;

typedef std::map<std::string, std::string> argtype;

/**
  Define a generic interface for derived classes to track interaction energy.
 */
class EnergyMap {
 public:
  //@{
  /** @name Arguments
    - default_value: set initial or cleared values to this (default: 0).
   */
  explicit EnergyMap(argtype args = argtype());
  explicit EnergyMap(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return the default value.
  double default_value() const { return default_value_; }

  /// Clear the interaction.
  virtual void clear(
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index);

  /// Update the interaction
  virtual double update(
      const double energy,
      const int part1_index,
      const int site1_index,
      const int site1_type,
      const int part2_index,
      const int site2_index,
      const int site2_type,
      const double squared_distance,
      const Position * pbc,
      const Configuration& config);

  /// Return true if the total interaction energy is stored (e.g., no criteria
  /// for inclusion.
  virtual bool is_queryable() const { return true; }

  /// Return the interaction energy of given site.
  virtual double energy(const int part1_index, const int site1_index) const;

  /// Precompute
  void precompute(Configuration * config);

  /// Revert any changes from perturbation of selection.
  virtual void revert(const Select& select) {}

  /// Finalize any changes from perturbation of selection.
  virtual void finalize(const Select& select) {}

  /// Return the total energy.
  virtual double total_energy() const;

  /**
    Add neighboring particles to selection which interact with node
    based on an energy less than or equal to the tolerance.
    The cluster also has positions taking into account periodic boundary
    conditions, which is why frame of reference is used recurisvely.
   */
  virtual void select_cluster(const NeighborCriteria& cluster_criteria,
                              const Configuration& config,
                              const int particle_node,
                              Select * cluster,
                              const Position& frame_of_reference) const;

  /// Compare old and new maps to see if cluster has changed.
  /// This is useful for detailed balance with rigid cluster moves.
  virtual bool is_cluster_changed(const NeighborCriteria& cluster_criteria,
    const Select& select,
    const Configuration& config) const;

  /// Return the NeighborCriteria.
  virtual const NeighborCriteria& neighbor_criteria() const;

  /// Return the neighbors of target particle and site that are of given
  /// site index.
  virtual void neighbors(
    const NeighborCriteria& neighbor_criteria,
    const Configuration& config,
    const int target_particle,
    const int target_site,
    const int given_site_index,
    /// Return of the neighbors
    Select * neighbors,
    /// If 1, use newly computed map.
    const int new_map = 0) const;

  virtual void check(const Configuration& config) const {}
  virtual void is_equal(const EnergyMap& map) const {}
  virtual void clear() { data_.clear(); }

  // virtual const std::vector<double>& map(const int part1, const int part2,
  //   const int site1, const int site2) const;

  // Synchronize with another object of the same type.
  // Typically used with prefetch.
  virtual void synchronize_(const EnergyMap& map, const Select& perturbed);
  const SynchronizeData& data() const { return data_; }

  // serialization
  virtual std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const {
    serialize_energy_map_(ostr); }
  explicit EnergyMap(std::istream& istr);
  virtual std::shared_ptr<EnergyMap> create(std::istream& istr) const;
  virtual std::shared_ptr<EnergyMap> create(argtype * args) const;
  std::map<std::string, std::shared_ptr<EnergyMap> >& deserialize_map();
  std::shared_ptr<EnergyMap> deserialize(std::istream& istr);
  std::shared_ptr<EnergyMap> factory(const std::string name, argtype * args);
  virtual ~EnergyMap() {}

  //@}
 protected:
  std::string class_name_ = "EnergyMap";
  SynchronizeData data_;
  void serialize_energy_map_(std::ostream& ostr) const;

  virtual void resize_(const int part1_index,
                       const int site1_index,
                       const int part2_index,
                       const int site2_index);

  virtual std::vector<double> * smap_(const int part1_index,
                                      const int site1_index,
                                      const int part2_index,
                                      const int site2_index);
  virtual std::vector<double> * smap_new_(const int part1_index,
                                          const int site1_index,
                                          const int part2_index,
                                          const int site2_index);

  virtual const std::vector<std::vector<std::vector<std::vector<std::vector<
    double> > > > >& map() const;
  int dimen() const { return dimen_; }
  int site_max() const { return site_max_; }

  // temporary and non-serialized
  std::vector<int> stored_neighbors_;

 private:
  double default_value_;
  // HWH optimization, could make variation in size but hard to initialize
  int site_max_;  // largest number of sites in a particle.
  int dimen_ = -1;  // record dimension of pbc wrap
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_ENERGY_MAP_H_
