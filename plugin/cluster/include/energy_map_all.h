
#ifndef FEASST_CLUSTER_ENERGY_MAP_ALL_H_
#define FEASST_CLUSTER_ENERGY_MAP_ALL_H_

#include <vector>
#include "system/include/energy_map.h"
#include "configuration/include/neighbor_criteria.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Map energies with an inefficient but simple data structure.
  All pairwise interactions between particles and sites are stored, even when
  the interaction is zero.

  This implementation stores two versions of the map, the current map and the
  new map.
  Updates from perturbations change only the new map.
  If the perturbation is accepted, the updates are finalized into the current map.
  Otherwise, the new map is synchronized to the old map.
 */
class EnergyMapAll : public EnergyMap {
 public:
  explicit EnergyMapAll(argtype args = argtype());
  explicit EnergyMapAll(argtype * args);
  double energy(const int part1_index, const int site1_index) const override;
  void revert(const Select& select) override;
  void finalize(const Select& select) override;
  void select_cluster(const NeighborCriteria& neighbor_criteria,
                      const Configuration& config,
                      const int particle_node,
                      Select * cluster,
                      const Position& frame_of_reference) const override;
  bool is_cluster_changed(const NeighborCriteria& neighbor_criteria,
    const Select& select,
    const Configuration& config) const override;
  void neighbors(
    const NeighborCriteria& neighbor_criteria,
    const Configuration& config,
    const int target_particle,
    const int target_site,
    const int given_site_index,
    Select * neighbors,
    const int new_map = 0) const override;
  void check(const Configuration& config) const override;
  void synchronize_(const EnergyMap& map, const Select& perturbed) override;
  void clear() override;

  // serialization
  std::string class_name() const override { return class_name_; }
  std::shared_ptr<EnergyMap> create(std::istream& istr) const override {
    return std::make_shared<EnergyMapAll>(istr); }
  std::shared_ptr<EnergyMap> create(argtype * args) const override {
    return std::make_shared<EnergyMapAll>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit EnergyMapAll(std::istream& istr);
  virtual ~EnergyMapAll() {}

 protected:
  void serialize_energy_map_all_(std::ostream& ostr) const;
  void resize_(const int part1, const int site1, const int part2, const int site2) override;
  std::vector<double> * smap_(const int part1_index,
                              const int site1_index,
                              const int part2_index,
                              const int site2_index) override {
    return &((*map_())[part1_index][part2_index][site1_index][site2_index]); }
  std::vector<double> * smap_new_(const int part1_index,
                                  const int site1_index,
                                  const int part2_index,
                                  const int site2_index) override {
    return &((*map_new_())[part1_index][part2_index][site1_index][site2_index]); }
  const std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > >& map() const override { return data_.dble_6D()[0]; }
  const std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > >& map_new() const { return data_.dble_6D()[1]; }

 private:
  std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > * map_();
  std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > * map_new_();
  int part_max_() { return static_cast<int>(map().size()); }
  bool is_cluster_(const NeighborCriteria& neighbor_criteria,
                   const std::vector<std::vector<std::vector<double> > >& smap,
                   const int particle_index0,
                   const int particle_index1,
                   const Configuration& config,
                   Position * frame = NULL) const;
};

inline std::shared_ptr<EnergyMapAll> MakeEnergyMapAll(
    const argtype& args = argtype()) {
  return std::make_shared<EnergyMapAll>(args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_ENERGY_MAP_ALL_H_
