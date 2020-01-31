
#ifndef FEASST_SYSTEM_ENERGY_MAP_ALL_H_
#define FEASST_SYSTEM_ENERGY_MAP_ALL_H_

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
  void clear(
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index) override;
  double update(
      const double energy,
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index,
      const double squared_distance,
      const Position * pbc) override;
  double query(
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index) override;
  void precompute(Configuration * config) override;
  void prep_for_revert(const Select& select) override;
  void revert(const Select& select) override;
  void remove_particles(const Select& select) override;
  double total_energy() const override;
  void select_cluster(const ClusterCriteria * cluster_criteria,
                      const Configuration& config,
                      const int particle_node,
                      SelectPosition * cluster,
                      const Position& frame_of_reference) const override;

  // serialization
  std::string class_name() const override { return class_name_; }
  std::shared_ptr<EnergyMap> create(std::istream& istr) const override {
    return std::make_shared<EnergyMapAll>(istr); }
  void serialize(std::ostream& ostr) const override;
  EnergyMapAll(std::istream& istr);
  virtual ~EnergyMapAll() {}

 private:
  const std::string class_name_ = "EnergyMapAll";
  std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > map_, map_old_;
  // HWH optimization, could make variation in size but hard to initialize
  int site_max_;  // largest number of sites in a particle.
  int dimen_ = -1; // record dimension of pbc wrap
  int part_max_() { return static_cast<int>(map_.size()); }
  void resize_(const int part1, const int site1, const int part2, const int site2);
  bool is_cluster_(const ClusterCriteria * cluster_criteria,
                   const std::vector<std::vector<std::vector<double> > >& smap,
                   Position * frame) const;
};

inline std::shared_ptr<EnergyMapAll> MakeEnergyMapAll(
    const argtype& args = argtype()) {
  return std::make_shared<EnergyMapAll>(args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_ENERGY_MAP_ALL_H_
