
#ifndef FEASST_SYSTEM_ENERGY_MAP_H_
#define FEASST_SYSTEM_ENERGY_MAP_H_

#include <string>
#include <vector>
#include <memory>
#include <map>
#include "utils/include/arguments.h"
#include "configuration/include/configuration.h"
#include "configuration/include/select.h"
#include "system/include/neighbor_criteria.h"

namespace feasst {

/**
  Define a generic interface for derived classes to track interaction energy.
 */
class EnergyMap {
 public:
  EnergyMap(
    /**
      args:
      - default_value: set initial or cleared values to this.
     */
    const argtype& args = argtype()) {
    Arguments args_(args);
    default_value_ = args_.key("default_value").dflt("0.").dble();
  }

  double default_value() const { return default_value_; }

  virtual void clear(
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index);
  virtual double update(
      const double energy,
      const int part1_index,
      const int site1_index,
      const int site1_type,
      const int part2_index,
      const int site2_index,
      const int site2_type,
      const double squared_distance,
      const Position * pbc);
  virtual bool is_queryable() const { return true; }
  double query(
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index);
  void precompute(Configuration * config);

  // HWH move to a finalize instead of revert-heavy stance
  // update partial map
  // Don't update full map until ...
  // finalize(selection)
  /* HWH
    For reverting, consider two different maps: total, and partial.
    The partial is zero except for the interactions of the selection.
    If move is accepted (finalize?) then total is replaced with partial.
    But theres no quick way to skip over the zeros (except using selection?)
    Or maybe its same efficiency to have old and new...
  */
  virtual void prep_for_revert(const Select& select) {}
  virtual void revert(const Select& select) {}
  virtual void finalize(const Select& select) {}
  double total_energy() const;

  /**
    Add neighboring particles to selection which interact with node
    based on an energy less than or equal to the tolerance.
    The cluster also has positions taking into account periodic boundary
    conditions, which is why frame of reference is used recurisvely.
   */
  virtual void select_cluster(const NeighborCriteria * cluster_criteria,
                              const Configuration& config,
                              const int particle_node,
                              Select * cluster,
                              const Position& frame_of_reference) const;

  /// Compare old and new maps to see if cluster has changed.
  /// This is useful for detailed balance with rigid cluster moves.
  virtual bool is_cluster_changed(const NeighborCriteria * cluster_criteria,
    const Select& select,
    const Configuration& config) const;

  /// Return the NeighborCriteria.
  virtual const NeighborCriteria * neighbor_criteria() const;

  /// Return a random neighboring site of target_site in target_particle.
  /// This interface was designed for use by AVB methods.
  virtual void neighbors(
    const NeighborCriteria * neighbor_criteria,
    const Configuration& config,
    const int target_particle,
    const int target_site,
    /// random_site is the given site index.
    const int random_site,
    Random * random,
    Select * neighbors) const;

  virtual void check() const {}

  // serialization
  virtual std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const {
    serialize_energy_map_(ostr); }
  explicit EnergyMap(std::istream& istr);
  virtual std::shared_ptr<EnergyMap> create(std::istream& istr) const = 0;
  std::map<std::string, std::shared_ptr<EnergyMap> >& deserialize_map();
  std::shared_ptr<EnergyMap> deserialize(std::istream& istr);
  virtual ~EnergyMap() {}

 protected:
  std::string class_name_ = "EnergyMap";
  void serialize_energy_map_(std::ostream& ostr) const;

  virtual void resize_(const int part1_index,
                       const int site1_index,
                       const int part2_index,
                       const int site2_index) = 0;

  virtual std::vector<double> * smap_(const int part1_index,
                                      const int site1_index,
                                      const int part2_index,
                                      const int site2_index) = 0;
  virtual std::vector<double> * smap_new_(const int part1_index,
                                          const int site1_index,
                                          const int part2_index,
                                          const int site2_index) = 0;

  virtual const std::vector<std::vector<std::vector<std::vector<std::vector<
    double> > > > >& map() const = 0;
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
