
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

  virtual void clear(
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index) = 0;
  virtual double update(
      const double energy,
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index,
      const double squared_distance,
      const Position * pbc) = 0;
  virtual double query(
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index) = 0;
  virtual void precompute(Configuration * config) {}

  /* HWH
    For reverting, consider two different maps: total, and partial.
    The partial is zero except for the interactions of the selection.
    If move is accepted (finalize?) then total is replaced with partial.
    But theres no quick way to skip over the zeros (except using selection?)
    Or maybe its same efficiency to have old and new...
  */
  virtual void prep_for_revert(const Select& select) {}
  virtual void revert(const Select& select) {}
  virtual void remove_particles(const Select& select) {}
  virtual double total_energy() const { FATAL("not implemented"); }

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

  // serialization
  virtual std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const {
    serialize_energy_map_(ostr); }
  EnergyMap(std::istream& istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(version == 945, "version mismatch: " << version);
    feasst_deserialize(&default_value_, istr);
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
  }

 private:
  double default_value_;
  const std::string class_name_ = "EnergyMap";
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_ENERGY_MAP_H_
