
#ifndef FEASST_SYSTEM_ENERGY_MAP_H_
#define FEASST_SYSTEM_ENERGY_MAP_H_

#include "utils/include/arguments.h"
#include "system/include/model.h"
#include "configuration/include/configuration.h"
#include "system/include/select_list.h"
#include "system/include/cluster_criteria.h"

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

  void clear(
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index);
  virtual double update(
      const double energy,
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index,
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
  virtual void select_cluster(const ClusterCriteria * cluster_criteria,
                              const Configuration& config,
                              const int particle_node,
                              SelectPosition * cluster,
                              const Position& frame_of_reference) const {
    FATAL("not implemented"); }

  /// Compare old and new maps to see if cluster has changed.
  /// This is useful for detailed balance with rigid cluster moves.
  virtual bool is_cluster_changed(const ClusterCriteria * cluster_criteria,
    const Select& select) const { FATAL("not implemented"); }

  virtual void check() const {}

  // serialization
  virtual std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const {
    serialize_energy_map_(ostr); }
  EnergyMap(std::istream& istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(version == 945, "version mismatch: " << version);
    feasst_deserialize(&default_value_, istr);
    feasst_deserialize(&site_max_, istr);
    feasst_deserialize(&dimen_, istr);
  }
  virtual std::shared_ptr<EnergyMap> create(std::istream& istr) const = 0;
  std::map<std::string, std::shared_ptr<EnergyMap> >& deserialize_map();
  std::shared_ptr<EnergyMap> deserialize(std::istream& istr) {
    return template_deserialize(deserialize_map(), istr); }
  virtual ~EnergyMap() {}

 protected:
  void serialize_energy_map_(std::ostream& ostr) const {
    ostr << class_name_ << " ";
    feasst_serialize_version(945, ostr);
    feasst_serialize(default_value_, ostr);
    feasst_serialize(site_max_, ostr);
    feasst_serialize(dimen_, ostr);
  }

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

  virtual const std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > >& map() const = 0;
  int dimen() const { return dimen_; }
  int site_max() const { return site_max_; }

 private:
  double default_value_;
  const std::string class_name_ = "EnergyMap";
  // HWH optimization, could make variation in size but hard to initialize
  int site_max_;  // largest number of sites in a particle.
  int dimen_ = -1; // record dimension of pbc wrap
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_ENERGY_MAP_H_
