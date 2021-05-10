#include <algorithm>
#include "utils/include/utils.h"  // find_in_list
#include "utils/include/io.h"
#include "utils/include/serialize.h"
#include "configuration/include/select.h"
#include "configuration/include/configuration.h"
#include "cluster/include/energy_map_neighbor.h"

namespace feasst {

class MapEnergyMapNeighbor {
 public:
  MapEnergyMapNeighbor() {
    EnergyMapNeighbor().deserialize_map()["EnergyMapNeighbor"] =
      MakeEnergyMapNeighbor();
  }
};

static MapEnergyMapNeighbor mapper_ = MapEnergyMapNeighbor();

EnergyMapNeighbor::EnergyMapNeighbor(const argtype& args) : EnergyMap(args) {
  class_name_ = "EnergyMapNeighbor";
}

EnergyMapNeighbor::EnergyMapNeighbor(std::istream& istr) : EnergyMap(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2947, "mismatch:" << version);
}

void EnergyMapNeighbor::serialize_energy_map_neighbor_(std::ostream& ostr) const {
  serialize_energy_map_(ostr);
  feasst_serialize_version(2947, ostr);
}

void EnergyMapNeighbor::serialize(std::ostream& ostr) const {
  serialize_energy_map_neighbor_(ostr);
}

std::vector<double> * EnergyMapNeighbor::smap_(const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index) {
  FATAL("not implemented");
}

std::vector<double> * EnergyMapNeighbor::smap_new_(const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index) {
  FATAL("not implemented");
}

double EnergyMapNeighbor::update(
    const double energy,
    const int part1_index,
    const int site1_index,
    const int site1_type,
    const int part2_index,
    const int site2_index,
    const int site2_type,
    const double squared_distance,
    const Position * pbc) {
  if (energy != 0.) {
    std::vector<double> * mn1 =
      find_or_add_(site2_index,
      find_or_add_(part2_index,
      find_or_add_(site1_index,
      find_or_add_(part1_index, map_new_()))));
    if (mn1->size() == 0) *mn1 = std::vector<double>(5, 0);
    (*mn1)[0] = energy;
    (*mn1)[1] = squared_distance;
    if (pbc->dimension() > 0) {
      for (int dim = 0; dim < dimen(); ++dim) {
        (*mn1)[2 + dim] = pbc->coord(dim);
      }
    }
  }
  finalizable_ = true;
  return energy;
}

void EnergyMapNeighbor::revert(const Select& select) {
  map_new_()->clear();
}

void EnergyMapNeighbor::sort_map_new_() {
  std::sort(map_new_()->begin(), map_new_()->end());
  for (std::pair<int, mn4type>& mn4 : *map_new_()) {
    std::sort(mn4.second.begin(), mn4.second.end());
    for (std::pair<int, map3type>& mn3 : mn4.second) {
      std::sort(mn3.second.begin(), mn3.second.end());
      for (std::pair<int, map2type>& mn2 : mn3.second) {
        std::sort(mn2.second.begin(), mn2.second.end());
      }
    }
  }
  DEBUG("sorted new map: " << map_new_str());
}

void EnergyMapNeighbor::size_map_() {
  for (const std::pair<int, mn4type>& mn4 : *map_new_()) {
    const int part1 = mn4.first;
    DEBUG("part1 " << part1);
    find_or_add_(part1, map_());
    DEBUG("mapsize " << map_()->size());
    for (const auto& mn3 : mn4.second) {
      const int site1 = mn3.first;
      DEBUG("site1 " << site1);
      find_or_add_(site1, &(*map_())[part1]);
      DEBUG("map4size " << (*map_())[part1].size());
    }
  }
  DEBUG("map size after init " << map_()->size())
}

void EnergyMapNeighbor::remove_from_map_nvt_(const Select& select) {
  // remove 1. part2 in selection not in new
  //        2. site2 in selection not in new
  for (int spindex = 0; spindex < select.num_particles(); ++spindex) {
    const int part1 = select.particle_index(spindex);
    DEBUG("part1 " << part1);
    if (part1 < static_cast<int>(map_()->size())) {
      map4type * map4 = &(*map_())[part1];
      mn4type * mn4 = find_or_add_(part1, map_new_());
      for (const int site1 : select.site_indices(spindex)) {
        if (site1 < static_cast<int>(map4->size())) {
          DEBUG("site1 " << site1);
          map3type * map3 = &(*map4)[site1];
          map3type * mn3 = find_or_add_(site1, mn4);

          // remove part2 in selection not in new map
          map3type diff3;
          std::set_difference(map3->begin(), map3->end(),
                              mn3->begin(), mn3->end(), back_inserter(diff3));
          for (const auto& mn2 : diff3) {
            const int part2 = mn2.first;
            // skip double counted
            DEBUG(select.trial_state());
            if (select.trial_state() == 1 || part2 > part1) {
              DEBUG("removing part1/2: " << part1 << "/" << part2);
              int fpindex = -1;
              const bool found = find_in_list(part2, *map3, &fpindex);
              if (found) {
                map3->erase(map3->begin() + fpindex,
                            map3->begin() + fpindex + 1);
              } else {
                FATAL(part2 << " not found");
              }

              // for removed part2, also remove perturbed indices
              for (const auto& mn1 : mn2.second) {
                const int site2 = mn1.first;
                map2type * map2inv = find_or_add_(part1,
                                     find_or_add_(site2,
                                     find_or_add_(part2, map_())));
                int fsindex = -1;
                const bool found = find_in_list(site1, *map2inv, &fsindex);
                if (found) {
                  map2inv->erase(map2inv->begin() + fsindex,
                                 map2inv->begin() + fsindex + 1);
                } else {
                  FATAL(part2 << " not found");
                }
              }
            }
          }

          // find part2 in selection and old map in order to compare sites
          map3type intersect3;
          std::set_intersection(map3->begin(), map3->end(),
                                mn3->begin(), mn3->end(),
                                back_inserter(intersect3));
          for (const auto& mn2 : intersect3) {
            const int part2 = mn2.first;
            if (select.trial_state() == 1 || part2 > part1) {
              // find site2 in selection not in new map
              map2type * map2 = find_or_add_(part2, map3);
              map2type diff2;
              std::set_difference(map2->begin(), map2->end(),
                                  mn2.second.begin(), mn2.second.end(),
                                  back_inserter(diff2));
              for (const auto& missing1 : diff2) {
                const int site2 = missing1.first;
                int fsindex = -1;
                bool found = find_in_list(site2, *map2, &fsindex);
                ASSERT(found, "err");
                map2->erase(map2->begin() + fsindex,
                            map2->begin() + fsindex + 1);
                // for removed site2, also remove perturbed index
                map2type * map2inv = find_or_add_(part1,
                                     find_or_add_(site2,
                                     find_or_add_(part2, map_())));
                fsindex = -1;
                found = find_in_list(site1, *map2inv, &fsindex);
                ASSERT(found, "err");
                map2inv->erase(map2inv->begin() + fsindex,
                               map2inv->begin() + fsindex + 1);
              }
            }
          }
        }
      }
    }
  }
}

void EnergyMapNeighbor::add_to_map_nvt_() {
  // add 1. part2 in new map not in old map
  //     2. site2 in new map not in old map
  DEBUG("map new size: " << map_new_()->size());
  for (std::pair<int, mn4type>& mn4 : *map_new_()) {
    const int part1 = mn4.first;
    DEBUG("part1 " << part1);
    map4type * map4 = find_or_add_(part1, map_());
    for (std::pair<int, map3type>& mn3 : mn4.second) {
      const int site1 = mn3.first;
      map4 = &(*map_())[part1];
      // resizing may make map4 pointer relocate, same with map3
      map3type * map3 = find_or_add_(site1, map4);

      // find part2 in new map not in old map
      map3type diff3;
      std::set_difference(mn3.second.begin(), mn3.second.end(),
                          map3->begin(), map3->end(), back_inserter(diff3));
      for (const auto& mn2 : diff3) {
        const int part2 = mn2.first;
        DEBUG("adding part1/2: " << part1 << "/" << part2);
        map3 = &(*map_())[part1][site1];
        map2type * map2 = find_or_add_(part2, map3);
        *map2 = mn2.second;
        std::sort(map3->begin(), map3->end());

        // for newly added part2, also add perturbed indices of each site
        for (const auto& mn1 : mn2.second) {
          const int site2 = mn1.first;
          map2type * map2inv = find_or_add_(part1,
                               find_or_add_(site2,
                               find_or_add_(part2, map_())));
          map1type * map1inv = find_or_add_(site1, map2inv);
          *map1inv = mn1.second;
          invert_pbcs_(map1inv);
          std::sort(map2inv->begin(), map2inv->end());
        }
      }

      // find part2 in new map and old map in order to compare the sites
      map3type intersect3;
      std::set_intersection(mn3.second.begin(), mn3.second.end(),
                            map3->begin(), map3->end(),
                            back_inserter(intersect3));
      for (const auto& mn2 : intersect3) {
        const int part2 = mn2.first;
        map2type * map2 = find_or_add_(part2, map3);
        // find site2 in new map not in old map
        map2type diff2;
        std::set_difference(mn2.second.begin(), mn2.second.end(),
                            map2->begin(), map2->end(), back_inserter(diff2));
        for (const auto& mn1 : diff2) {
          const int site2 = mn1.first;
          DEBUG("adding site2: " << site2);
          map1type * map1 = find_or_add_(site2, map2);
          *map1 = mn1.second;
          std::sort(map2->begin(), map2->end());

          // for newly added site2, also add perturbed index
          // HWH copied from above
          map2type * map2inv = find_or_add_(part1,
                               find_or_add_(site2,
                               find_or_add_(part2, map_())));
          map1type * map1inv = find_or_add_(site1, map2inv);
          *map1inv = mn1.second;
          invert_pbcs_(map1inv);
          std::sort(map2inv->begin(), map2inv->end());
        }
      }
    }
  }

  DEBUG("map after adding before removing: " << map_str());
  DEBUG("new map after adding before removing: " << map_new_str());
}

void EnergyMapNeighbor::remove_particle_from_map_(const Select& select) {
  //ASSERT(select.num_particles() == 1, "multi-particle not implemented");
  for (int spindex = 0; spindex < select.num_particles(); ++spindex) {
    const int part1 = select.particle_index(spindex);
    DEBUG("part1 " << part1);
    if (part1 < static_cast<int>(map_()->size())) {
      map4type * map4 = &(*map_())[part1];
      for (const int site1 : select.site_indices(spindex)) {
        if (site1 < static_cast<int>(map4->size())) {
          DEBUG("site1 " << site1);
          map3type * map3 = &(*map4)[site1];
          for (const std::pair<int, map2type>& map2 : *map3) {
            const int part2 = map2.first;
            DEBUG("part2 " << part2);
            for (const std::pair<int, map1type>& map1 : map2.second) {
              const int site2 = map1.first;
              DEBUG("site2 " << site2);
              map2type * map2inv = find_or_add_(part1,
                                   find_or_add_(site2,
                                   find_or_add_(part2, map_())));
              int fsindex = -1;
              const bool found = find_in_list(site2, *map2inv, &fsindex);
              if (found) {
                map2inv->erase(map2inv->begin() + fsindex,
                               map2inv->begin() + fsindex + 1);
              } else {
                FATAL(site2 << " not found");
              }
            }
          }
        }
      }
      (*map_())[part1] = map4type(); // erase
    }
  }
}

void EnergyMapNeighbor::add_particle_to_map_() {
  // Add "new" to map
  DEBUG("map new size: " << map_new_()->size());
  for (const std::pair<int, mn4type>& mn4 : *map_new_()) {
    const int part1 = mn4.first;
    DEBUG("part1 " << part1);
    map4type * map4 = find_or_add_(part1, map_());
    for (const std::pair<int, map3type>& mn3 : mn4.second) {
      const int site1 = mn3.first;
      DEBUG("site1 " << site1);
      map3type * map3 = find_or_add_(site1, map4);
      for (const std::pair<int, map2type>& mn2 : mn3.second) {
        const int part2 = mn2.first;
        DEBUG("part2 " << part2);
        map2type * map2 = find_or_add_(part2, map3);
        for (const std::pair<int, map1type>& mn1 : mn2.second) {
          const int site2 = mn1.first;
          DEBUG("site2 " << site2);
          map1type * map1 = find_or_add_(site1, map2);
          *map1 = mn1.second;

          // HWH copied from above.. make function
          // perturb indices and add
          map2type * map2inv = find_or_add_(part1,
                               find_or_add_(site2,
                               find_or_add_(part2, map_())));
          map1type * map1inv = find_or_add_(site1, map2inv);
          *map1inv = mn1.second;
          invert_pbcs_(map1inv);
          std::sort(map2inv->begin(), map2inv->end());
        }
      }
    }
  }
}

void EnergyMapNeighbor::finalize(const Select& select) {
  if (!finalizable_) return;
  DEBUG("map: " << map_str());
  DEBUG("new map: " << map_new_str());
  DEBUG("perturbed: " << select.str());

  DEBUG("map size " << map_()->size())
  DEBUG("new map size " << map_new_()->size())

  sort_map_new_();
  size_map_();

  if (select.trial_state() == -1 || select.trial_state() == 1) {
    remove_from_map_nvt_(select);
  }

  if (select.trial_state() == -1 || select.trial_state() == 1) {
    add_to_map_nvt_();
  } else if (select.trial_state() == 2) {
    remove_particle_from_map_(select);
  } else if (select.trial_state() == 3) {
    add_particle_to_map_();
  } else {
    FATAL("unrecognized trial state: " << select.trial_state());
  }
  map_new_()->clear();
}

void EnergyMapNeighbor::select_cluster(
    const NeighborCriteria& neighbor_criteria,
    const Configuration& config,
    const int particle_node,
    Select * cluster,
    const Position& frame_of_reference) const {
  DEBUG("particle_node " << particle_node);
  DEBUG("map size " << const_map_().size());
  if (particle_node < static_cast<int>(const_map_().size())) {
    const map4type& map4 = const_map_()[particle_node];
    for (int site1 = 0; site1 < static_cast<int>(map4.size()); ++site1) {
      const map3type& map3 = map4[site1];
      for (const std::pair<int, map2type>& map2 : map3) {
        const int part2_index = map2.first;
        /// This part is copy-pased from EnergyMapAll, except for is_cluster_
        // if part2 isn't already in the cluster
        // and part2 satistifies cluster criteria,
        // then recurively add part2 as a new node
        if (!find_in_list(part2_index, cluster->particle_indices())) {
          Position frame;
          if (is_cluster_(neighbor_criteria,
                          particle_node,
                          site1,
                          part2_index,
                          config,
                          true,  // old map
                          &frame)) {
            DEBUG("FOR " << frame_of_reference.str());
            frame.add(frame_of_reference);
            const Particle& part = config.select_particle(part2_index);
            cluster->add_particle(part, part2_index);
            cluster->load_positions_of_last(part, frame);
            DEBUG("frame: " << frame.str());
            DEBUG("added: " << cluster->site_positions().back().back().str());
            select_cluster(neighbor_criteria, config, part2_index, cluster, frame);
          }
        }
      }
    }
  }
}

//const map3type&
typedef std::vector<double> map1type;
typedef std::vector<std::pair<int, map1type> > map2type;
typedef std::vector<std::pair<int, map2type> > map3type;
const map3type * EnergyMapNeighbor::find_map3_(const int particle_index1,
    const int site_index1,
    const bool old) const {
  if (old) {
    if (particle_index1 < static_cast<int>(const_map_().size())) {
      const map4type& map4 = const_map_()[particle_index1];
      if (site_index1 < static_cast<int>(map4.size())) {
        return const_cast<const map3type *>(&map4[site_index1]);
      }
    }
  } else {
    int findex = -1;
    if (find_in_list(particle_index1, const_map_new_(), &findex)) {
      const mn4type& map4 = const_map_new_()[findex].second;
      if (find_in_list(site_index1, map4, &findex)) {
        return const_cast<const map3type *>(&map4[findex].second);
      }
    }
  }
  return NULL;
}

const map2type * EnergyMapNeighbor::find_map2_(const int part1,
    const int site1,
    const int part2,
    const bool old) const {
  const map3type * map3 = find_map3_(part1, site1, old);
  int findex = -1;
  if (find_in_list(part2, *map3, &findex)) {
    return const_cast<const map2type *>(&(*map3)[findex].second);
  }
  return NULL;
}

bool EnergyMapNeighbor::is_cluster_(const NeighborCriteria& neighbor_criteria,
    const int particle_index1,
    const int site_index1,
    const int particle_index2,
    const Configuration& config,
    const bool old,
    Position * frame) const {
  const map2type * map2 = find_map2_(particle_index1, site_index1,
                                     particle_index2, old);
  if (map2) {
    const int site_type1 =
      config.select_particle(particle_index1).site(site_index1).type();
    for (const std::pair<int, map1type>& map1 : *map2) {
      const int site2 = map1.first;
      const int site_type2 =
        config.select_particle(particle_index2).site(site2).type();

      // below was copied from EnergyMapAll, but with map1.second
      if (neighbor_criteria.is_accepted(map1.second[0], map1.second[1],
                                         site_type1, site_type2)) {
        if (frame) {
          frame->set_to_origin(dimen());
          for (int dim = 0; dim < dimen(); ++dim) {
            frame->set_coord(dim, -1.*map1.second[2 + dim]);
          }
        }
        return true;
      }
    }
  }
  return false;
}

bool EnergyMapNeighbor::is_cluster_changed(
    const NeighborCriteria& neighbor_criteria,
    const Select& select,
    const Configuration& config) const {
  for (int spindex = 0; spindex < select.num_particles(); ++spindex) {
    const int p1 = select.particle_index(spindex);
    for (const int s1 : select.site_indices(spindex)) {
      const map3type * newmap = find_map3_(p1, s1, false);
      const map3type * oldmap = find_map3_(p1, s1, true);
      if (!newmap && !oldmap) {
        return true;
      } else if (newmap && oldmap) {
        if (*newmap != *oldmap) {
          return true;
        }
      }
    }
  }
  return false;
}

void EnergyMapNeighbor::resize_(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index) { FATAL("not impl"); }

double EnergyMapNeighbor::total_energy() const {
  double en = 0.;
  for (const map4type& map4 : const_map_()) {
    for (const auto& map3 : map4) {
      for (const auto& map2 : map3) {
        for (const auto& map1 : map2.second) {
          en += map1.second[0];//map3.second[0][0][0];
        }
      }
    }
  }
  return 0.5*en;
}

void EnergyMapNeighbor::check() const {
  // see if all neighbors are double counted
  for (int part1 = 0; part1 < static_cast<int>(const_map_().size()); ++part1) {
    const map4type& m4 = const_map_()[part1];
    for (int site1 = 0; site1 < static_cast<int>(m4.size()); ++site1) {
      const map3type& m3 = m4[site1];
      for (int pneigh = 0; pneigh < static_cast<int>(m3.size()); ++pneigh) {
        const int part2 = m3[pneigh].first;
        const map2type& m2 = m3[pneigh].second;
        if (m2.size() > 0) {
          for (int sneigh = 0; sneigh < static_cast<int>(m2.size()); ++sneigh) {
            const int site2 = m2[sneigh].first;
            const map1type& m1 = m2[sneigh].second;
            if (m1.size() > 0) {
              int tmp;
              bool found = find_in_list(part1, const_map_()[part2][site2], &tmp);
              ASSERT(found, "unmatched pair part1: " << part1 << " part2: " << part2
                << " map: " << map_str());
              found = find_in_list(site1, const_map_()[part2][site2][tmp].second, &tmp);
              ASSERT(found, "unmatched pair part1: " << part1 << " site1: " <<
                site1 << " part2 " << part2 << " site2: " << site2
                << " map: " << map_str());
            }
          }
        }
      }
    }
  }
}

void EnergyMapNeighbor::neighbors(
    const NeighborCriteria& neighbor_criteria,
    const Configuration& config,
    const int target_particle,
    const int target_site,
    const int given_site_index,
    Select * neighbors,
    const int new_map) const {
  neighbors->clear();
  const Site& site0 = config.select_particle(target_particle).site(target_site);
  const int site_type0 = site0.type();
  DEBUG("target_particle " << target_particle);
  DEBUG("target_site " << target_site);
  //const vec4 * map4 = const_cast<vec4 * const>(&map()[target_particle]);
  const map3type * map3 = NULL;
  ASSERT(!map3, "er");
  map3 = find_map3_(target_particle, target_site, !new_map);
  if (map3) {
    DEBUG(map_str(*map3));
    for (int pindex = 0; pindex < static_cast<int>(map3->size()); ++pindex) {
      const int part2 = (*map3)[pindex].first;
      DEBUG("part2 " << part2);
      const map2type& map2 = (*map3)[pindex].second;
      int findex = -1;
      if (find_in_list(given_site_index, map2, &findex)) {
        const Site& site1 = config.select_particle(part2).site(given_site_index);
        const int site_type1 = site1.type();
        const map1type& map1 = map2[findex].second;
        DEBUG("site2: " << map2[findex].first);
        if (neighbor_criteria.is_accepted(map1[0], map1[1],
                                          site_type0, site_type1)) {
          neighbors->add_site(part2, given_site_index);
        }
      }
    }
  }
}

std::string EnergyMapNeighbor::map_str(const map3type& map3) const {
  std::stringstream ss;
  for (const std::pair<int, map2type>& map2 : map3) {
    ss << "p2 " << map2.first << ":{";
    for (const std::pair<int, map1type>& map1 : map2.second) {
      ss << "s2 " << map1.first << ":(" << feasst_str(map1.second) << "),";
    }
    ss << "},";
  }
  return ss.str();
}

std::string EnergyMapNeighbor::map_new_str() const {
  std::stringstream ss;
  ss << std::endl;
  for (const std::pair<int, mn4type>& map4: const_map_new_()) {
    ss << "p1 " << map4.first << ":|";
    for (const std::pair<int, map3type>& map3 : map4.second) {
      ss << "s1 " << map3.first << ":[" << map_str(map3.second) << "],";
    }
    ss << "|," << std::endl;
  }
  return ss.str();
}

std::string EnergyMapNeighbor::map_str() const {
  std::stringstream ss;
  ss << std::endl;
  for (int part1 = 0; part1 < static_cast<int>(const_map_().size()); ++part1) {
    ss << "p1 " << part1 << ":|";
    const map4type& map4 = const_map_()[part1];
    for (int site1 = 0; site1 < static_cast<int>(map4.size()); ++site1) {
      const map3type& map3 = map4[site1];
      ss << "s1 " << site1 << ":[" << map_str(map3) << "],";
    }
    ss << "|," << std::endl;
  }
  return ss.str();
}

const std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > >& EnergyMapNeighbor::map() const { FATAL("not impl"); }

typedef std::vector<std::pair<int, map3type> > mn4type;
typedef std::vector<map3type> map4type;
std::vector<map4type> * EnergyMapNeighbor::map_() {
  return data_.get_vvvpvpv();
}

const std::vector<map4type>& EnergyMapNeighbor::const_map_() const {
  return const_cast<const std::vector<map4type>&>(data_.get_const_vvvpvpv());
}

std::vector<std::pair<int, mn4type> > * EnergyMapNeighbor::map_new_() {
  return data_.get_vpvpvpvpv();
}

const std::vector<std::pair<int, mn4type> >& EnergyMapNeighbor::const_map_new_() const {
  return const_cast<std::vector<std::pair<int, mn4type> >&>(data_.get_const_vpvpvpvpv());
}

}  // namespace feasst
