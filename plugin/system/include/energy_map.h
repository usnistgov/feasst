
#ifndef FEASST_SYSTEM_ENERGY_MAP_H_
#define FEASST_SYSTEM_ENERGY_MAP_H_

#include <vector>
#include "system/include/model.h"
#include "configuration/include/configuration.h"

namespace feasst {

/**
  Define a generic interface for derived classes to track interaction energy.
 */
class EnergyMap {
 public:
  EnergyMap() {}
  virtual double update(
      const double energy,
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index) {
    return energy;
  }
  virtual void precompute(Configuration * config) {}

  /* HWH
    For reverting, consider two different maps: total, and partial.
    The partial is zero except for the interactions of the selection.
    If move is accepted (finalize?) then total is replaced with partial.
    But theres no quick way to skip over the zeros (except using selection?)
    Or maybe its same efficiency to have old and new...
  */
  virtual void prep_for_revert() {}
  virtual void revert() {}

  virtual const std::vector<std::vector<std::vector<std::vector<double> > > >&
    map() const { ERROR("not implemented"); }

  virtual double total() const { ERROR("not implemented"); }

  // serialization
  virtual std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const {
    serialize_energy_map_(ostr); }
  EnergyMap(std::istream& istr) { feasst_deserialize_version(istr); }
  virtual std::shared_ptr<EnergyMap> create(std::istream& istr) const {
    return std::make_shared<EnergyMap>(istr); }
  std::map<std::string, std::shared_ptr<EnergyMap> >& deserialize_map();
  std::shared_ptr<EnergyMap> deserialize(std::istream& istr) {
    return template_deserialize(deserialize_map(), istr); }
  virtual ~EnergyMap() {}

 protected:
  void serialize_energy_map_(std::ostream& ostr) const {
    ostr << class_name_ << " ";
    feasst_serialize_version(1, ostr);
  }

 private:
  const std::string class_name_ = "EnergyMap";
};

class EnergyMapAll : public EnergyMap {
 public:
  EnergyMapAll() {}
  double update(
      const double energy,
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index) override;
  void precompute(Configuration * config) override;
  void prep_for_revert() override { map_old_ = map_; }
  void revert() override { map_ = map_old_; }
  const std::vector<std::vector<std::vector<std::vector<double> > > >& map()
    const override { return map_; }
  double total() const override { return 0.5*sum(map()); }

  // serialization
  std::string class_name() const override { return class_name_; }
  std::shared_ptr<EnergyMap> create(std::istream& istr) const override {
    return std::make_shared<EnergyMapAll>(istr); }
  void serialize(std::ostream& ostr) const override;
  EnergyMapAll(std::istream& istr);
  virtual ~EnergyMapAll() {}

 private:
  const std::string class_name_ = "EnergyMapAll";
  std::vector<std::vector<std::vector<std::vector<double> > > > map_, map_old_;
  // HWH optimization, could make variation in size but hard to initialize
  int site_max_;  // largest number of sites in a particle.
  int part_max_() { return static_cast<int>(map_.size()); }
};

inline std::shared_ptr<EnergyMapAll> MakeEnergyMapAll() {
  return std::make_shared<EnergyMapAll>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_ENERGY_MAP_H_
