
#ifndef FEASST_CLUSTER_ENERGY_MAP_ALL_H_
#define FEASST_CLUSTER_ENERGY_MAP_ALL_H_

#include <vector>
#include "utils/include/arguments.h"
#include "system/include/energy_map.h"
#include "system/include/cluster_criteria.h"

namespace feasst {

/**
  Map energies with an inefficient but simple data structure.
  All pairwise interactions between particles and sites are stored.
 */
class EnergyMapAll : public EnergyMap {
 public:
  EnergyMapAll(const argtype& args = argtype()) : EnergyMap(args) {}
  void revert(const Select& select) override;
  void finalize(const Select& select) override;
  void select_cluster(const ClusterCriteria * cluster_criteria,
                      const Configuration& config,
                      const int particle_node,
                      SelectPosition * cluster,
                      const Position& frame_of_reference) const override;
  bool is_cluster_changed(const ClusterCriteria * cluster_criteria,
    const Select& select) const override;
  void check() const override;

  // serialization
  std::string class_name() const override { return class_name_; }
  std::shared_ptr<EnergyMap> create(std::istream& istr) const override {
    return std::make_shared<EnergyMapAll>(istr); }
  void serialize(std::ostream& ostr) const override;
  EnergyMapAll(std::istream& istr);
  virtual ~EnergyMapAll() {}

 protected:
  void serialize_energy_map_all_(std::ostream& ostr) const;
  void resize_(const int part1, const int site1, const int part2, const int site2) override;
  std::vector<double> * smap_(const int part1_index,
                              const int site1_index,
                              const int part2_index,
                              const int site2_index) override {
    return &map_[part1_index][part2_index][site1_index][site2_index]; }
  std::vector<double> * smap_new_(const int part1_index,
                                  const int site1_index,
                                  const int part2_index,
                                  const int site2_index) override {
    return &map_new_[part1_index][part2_index][site1_index][site2_index]; }
  const std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > >& map() const override { return map_; }

 private:
  const std::string class_name_ = "EnergyMapAll";
  std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > map_, map_new_;
  int part_max_() { return static_cast<int>(map_.size()); }
  bool is_cluster_(const ClusterCriteria * cluster_criteria,
                   const std::vector<std::vector<std::vector<double> > > * smap,
                   Position * frame = NULL) const;
};

inline std::shared_ptr<EnergyMapAll> MakeEnergyMapAll(
    const argtype& args = argtype()) {
  return std::make_shared<EnergyMapAll>(args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_ENERGY_MAP_ALL_H_
