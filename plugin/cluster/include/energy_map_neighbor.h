
#ifndef FEASST_CLUSTER_ENERGY_MAP_NEIGHBOR_H_
#define FEASST_CLUSTER_ENERGY_MAP_NEIGHBOR_H_

#include <vector>
//#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "system/include/energy_map.h"
#include "configuration/include/neighbor_criteria.h"

namespace feasst {

typedef std::vector<double> map1type;
typedef std::vector<std::pair<int, map1type> > map2type;
typedef std::vector<std::pair<int, map2type> > map3type;
typedef std::vector<std::pair<int, map3type> > mn4type;
typedef std::vector<map3type> map4type;

/**
  Map only between particles that interact (e.g., non zero energy).
  This data structure is intended to better scale with more particles than
  EnergyMapAll.
  Although for small system sizes or large cutoffs, EnergyMapAll may be faster
  because it does not require sorting.

  This implementation also stores two versions of the map.
  The first version is the current map.
  The second version is an empty new map.
  Upon perturbation, the empty new map is populated with interactions.
  When a perturbation is accepted, finalize replaces interactions in current
  map with those in the new map.
  The new map is then emptied.
  When a perturbation is rejected, revert empties the new map.
 */
class EnergyMapNeighbor : public EnergyMap {
 public:
  explicit EnergyMapNeighbor(argtype args = argtype());
  explicit EnergyMapNeighbor(argtype * args);
  double energy(const int part1_index, const int site1_index) const override;
  double update(
      const double energy,
      const int part1_index,
      const int site1_index,
      const int site1_type,
      const int part2_index,
      const int site2_index,
      const int site2_type,
      const double squared_distance,
      const Position * pbc,
      const Configuration& config) override;
  void revert(const Select& select) override;
  void finalize(const Select& select) override;
  double total_energy() const override;
  void check(const Configuration& config) const override;
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

  /// Clear interaction does nothing, because new map begins empty.
  void clear(
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index) override {}

  const std::vector<map4type>& const_map_() const;
  const std::vector<std::pair<int, mn4type> >& const_map_new_() const;

  // serialization
  std::string class_name() const override { return class_name_; }
  std::shared_ptr<EnergyMap> create(std::istream& istr) const override {
    return std::make_shared<EnergyMapNeighbor>(istr); }
  std::shared_ptr<EnergyMap> create(argtype * args) const override {
    return std::make_shared<EnergyMapNeighbor>(args); }
  void serialize(std::ostream& ostr) const override;
  EnergyMapNeighbor(std::istream& istr);
  virtual ~EnergyMapNeighbor() {}

 protected:
  void serialize_energy_map_neighbor_(std::ostream& ostr) const;
  void resize_(const int part1, const int site1, const int part2, const int site2) override;
  std::vector<double> * smap_(const int part1_index,
                              const int site1_index,
                              const int part2_index,
                              const int site2_index) override;
  std::vector<double> * smap_new_(const int part1_index,
                                  const int site1_index,
                                  const int part2_index,
                                  const int site2_index) override;
  const std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > >& map() const override;

 private:
  /// map_[part1][site1][pneigh].first -> part2
  ///                           .second[sneigh1].first-> site1
  ///                                           .second-> en, rsq, pbcs
  std::vector<map4type> * map_();

  /// map_[pneigh1].first -> part1
  ///              .second[sneigh1].first -> site1
  ///                              .second[pneigh2].first -> part2
  ///                                              .second[sneigh1].first -> site2
  ///                                                              .second -> en, rsq, pbcs
  std::vector<std::pair<int, mn4type> > * map_new_();
  //std::vector<std::vector<std::vector<std::vector<std::pair<int, std::vector<double> > > > > > map_new_;

//  /// The first index is the particle index, which mirrors config
//  /// The second index is the list of particles that are neighbors.
//  std::vector<std::vector<int> > neighbor_;
//
//  /// As opposed to EnergyMapAll, the second index (for particles) corresponds
//  /// with the second index of neighbor_ above, instead of listing all particles
//  /// in the configuration.
//  std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > map_;
//
//  /// Temporary and not to be serialized or synchronized.
//  /// These two data structures are to be populated by update and emptied by
//  /// revert/finalize
//  /// The first index is for the list of newly updated particle indices.
//  std::vector<int> updated_;
//  std::vector<std::vector<int> > neighbor_new_;
//
//  std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > map_new_;

  int part_max_() const { return static_cast<int>(const_map_().size()); }
  bool is_cluster_(const NeighborCriteria& neighbor_criteria,
                   const int particle_index1,
                   const int site_index1,
                   const int particle_index2,
                   const Configuration& config,
                   const bool old = true,
                   Position * frame = NULL) const;

  template<class T>
  T * find_or_add_(const int sindex, std::vector<T> * list) {
    const int missing = sindex - static_cast<int>(list->size()) + 1;
    //DEBUG("missing: " << missing);
    for (int index = 0; index < missing; ++index) list->push_back(T());
    return &(*list)[sindex];
  }

//  template<class T>
//  const T& find_(const int sindex, const std::vector<T>& list) const {
//    ASSERT(sindex < static_cast<int>(list.size()), "sindex: " << sindex
//      << " size: " << list.size());
//    return list[sindex];
//  }

  template<class T>
  T * find_or_add_(const int sindex, std::vector<std::pair<int, T> > * list) {
    int findex = -1;
    //DEBUG("sindex " << sindex);
    //DEBUG("list size " << list->size());
    if (!find_in_list(sindex, *list, &findex)) {
      findex = static_cast<int>(list->size());
      list->push_back(std::pair<int, T>());
      (*list)[findex].first = sindex;
      (*list)[findex].second = T();
    }
    //DEBUG("findex " << findex);
    return &(*list)[findex].second;
  }

//  template<class T>
//  const T& find_(const int sindex,
//      const std::vector<std::pair<int, T> >& list) const {
//    int findex = -1;
//    const bool found = find_in_list(sindex, list, &findex);
//    if (!found) FATAL("sindex:" << sindex << " not found.");
//    DEBUG("sindex " << sindex);
//    DEBUG("list size " << list.size());
//    return list[findex].second;
//  }

  const map3type * find_map3_(const int part1,
    const int site1,
    const bool old = true) const;

  const map2type * find_map2_(const int part1,
    const int site1,
    const int part2,
    const bool old = true) const;

  // invert pbcs
  void invert_pbcs_(map1type * map1) {
    for (int index = 2; index < static_cast<int>(map1->size()); ++index) {
      (*map1)[index] *= -1;
    }
  }

  // DEBUG util
  std::string map_new_str() const;
  std::string map_str() const;
  std::string map_str(const map3type& map3) const;

  void sort_map_new_();
  void size_map_();
  void remove_from_map_nvt_(const Select& select);
  void add_to_map_nvt_();
  void remove_particle_from_map_(const Select& select);
  void add_particle_to_map_();

  // temporary and not serialized
  bool finalizable_ = false;
  std::vector<std::vector<double> > energy_;
};

inline std::shared_ptr<EnergyMapNeighbor> MakeEnergyMapNeighbor(
    const argtype& args = argtype()) {
  return std::make_shared<EnergyMapNeighbor>(args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_ENERGY_MAP_NEIGHBOR_H_
