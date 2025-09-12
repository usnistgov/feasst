#include "utils/include/utils.h"  // find_in_list
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/select.h"
#include "configuration/include/configuration.h"
#include "cluster/include/energy_map_all.h"

namespace feasst {

FEASST_MAPPER(EnergyMapAll,);

EnergyMapAll::EnergyMapAll(argtype * args) : EnergyMap(args) {
  class_name_ = "EnergyMapAll";
  clear(); // resize
}
EnergyMapAll::EnergyMapAll(argtype args) : EnergyMapAll(&args) {
  feasst_check_all_used(args);
}

void EnergyMapAll::clear() {
  EnergyMap::clear();
  data_.get_dble_6D()->resize(2);
}

EnergyMapAll::EnergyMapAll(std::istream& istr) : EnergyMap(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2810, "mismatch:" << version);
}

void EnergyMapAll::serialize_energy_map_all_(std::ostream& ostr) const {
  serialize_energy_map_(ostr);
  feasst_serialize_version(2810, ostr);
}

void EnergyMapAll::serialize(std::ostream& ostr) const {
  serialize_energy_map_all_(ostr);
}

void EnergyMapAll::resize_(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index) {
  const int part_max_old = static_cast<int>(map().size());
  if (part1_index >= part_max_old || part2_index >= part_max_old) {
    const int part_max = std::max(part1_index, part2_index) + 1;
    map_()->resize(part_max);
    map_new_()->resize(part_max);
    for (int i = 0; i < part_max; ++i) {
      (*map_())[i].resize(part_max);
      (*map_new_())[i].resize(part_max);
      for (int j = 0; j < part_max; ++j) {
        ASSERT(site_max() != 0, "wasn't precomputed");
        (*map_())[i][j].resize(site_max());
        (*map_new_())[i][j].resize(site_max());
        for (int k = 0; k < site_max(); ++k) {
          (*map_())[i][j][k].resize(site_max());
          (*map_new_())[i][j][k].resize(site_max());
          if (i >= part_max_old || j >= part_max_old) {
            for (int l = 0; l < site_max(); ++l) {
              ASSERT(dimen() != -1, "wasnt precomputed");
              (*map_())[i][j][k][l] = std::vector<double>(2 + dimen(), default_value());
              (*map_new_())[i][j][k][l] = map()[i][j][k][l];
            }
          }
        }
      }
    }
  }
}

void EnergyMapAll::revert(const Select& select) {
  if (select.num_particles() > 0) {
    const int pmax = select.particle_indices().back();
    const int smax = select.site_indices().back().back();
    resize_(pmax, smax, pmax, smax);
    for (int sel_index = 0; sel_index < select.num_particles(); ++sel_index) {
      const int p1 = select.particle_index(sel_index);
      for (const int s1 : select.site_indices(sel_index)) {
        for (int p2 = 0; p2 < static_cast<int>(map()[p1].size()); ++p2) {
          for (int s2 = 0; s2 < static_cast<int>(map()[p1][p2][s1].size()); ++s2) {
            if (select.trial_state() == 3) {
              // revert addition
              (*map_())[p1][p2][s1][s2] = std::vector<double>(2 + dimen(), default_value());
              (*map_())[p2][p1][s2][s1] = std::vector<double>(2 + dimen(), default_value());
              (*map_new_())[p1][p2][s1][s2] = std::vector<double>(2 + dimen(), default_value());
              (*map_new_())[p2][p1][s2][s1] = std::vector<double>(2 + dimen(), default_value());
            } else {
              (*map_new_())[p1][p2][s1][s2] = map()[p1][p2][s1][s2];
              (*map_new_())[p2][p1][s2][s1] = map()[p2][p1][s2][s1];
            }
          }
        }
      }
    }
  }
}

void EnergyMapAll::finalize(const Select& select) {
  //INFO("finalizing: " << select.str());
  if (select.num_particles() > 0) {
    const int pmax = select.particle_indices().back();
    const int smax = select.site_indices().back().back();
    resize_(pmax, smax, pmax, smax);
    for (int sel_index = 0; sel_index < select.num_particles(); ++sel_index) {
      const int p1 = select.particle_index(sel_index);
      for (const int s1 : select.site_indices(sel_index)) {
        for (int p2 = 0; p2 < static_cast<int>(map()[p1].size()); ++p2) {
          for (int s2 = 0; s2 < static_cast<int>(map()[p1][p2][s1].size()); ++s2) {
            TRACE("p1s1p2s2 " << p1 << " " << s1 << " " << p2 << " " << s2);
            if (select.trial_state() == 2) {
              // finalize removal
              (*map_())[p1][p2][s1][s2] = std::vector<double>(2 + dimen(), default_value());
              (*map_())[p2][p1][s2][s1] = std::vector<double>(2 + dimen(), default_value());
              (*map_new_())[p1][p2][s1][s2] = std::vector<double>(2 + dimen(), default_value());
              (*map_new_())[p2][p1][s2][s1] = std::vector<double>(2 + dimen(), default_value());
            } else {
              (*map_())[p1][p2][s1][s2] = map_new()[p1][p2][s1][s2];
              (*map_())[p2][p1][s2][s1] = map_new()[p2][p1][s2][s1];
            }
          }
        }
      }
    }
  }
}

void EnergyMapAll::select_cluster(const NeighborCriteria& neighbor_criteria,
                                  const Configuration& config,
                                  const int particle_node,
                                  Select * cluster,
                                  const Position& frame_of_reference) const {
  DEBUG("particle_node " << particle_node);
  DEBUG("map size " << map().size());
  for (int part2_index = 0;
       part2_index < static_cast<int>(map()[particle_node].size());
       ++part2_index) {
    DEBUG("part2_index " << part2_index);
    // if part2 isn't already in the cluster
    // and part2 satistifies cluster criteria,
    // then recurively add part2 as a new node
    if (!find_in_list(part2_index, cluster->particle_indices())) {
      Position frame;
      if (is_cluster_(neighbor_criteria,
                      map()[particle_node][part2_index],
                      particle_node,
                      part2_index,
                      config,
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

bool EnergyMapAll::is_cluster_(
    const NeighborCriteria& neighbor_criteria,
    const std::vector<std::vector<std::vector<double> > >& smap,
    const int particle_index0,
    const int particle_index1,
    const Configuration& config,
    Position * frame) const {
  for (int s0i = 0; s0i < static_cast<int>(smap.size()); ++s0i) {
    const Site& site0 = config.select_particle(particle_index0).site(s0i);
    const int site_type0 = site0.type();
    const std::vector<std::vector<double> >& map2 = smap[s0i];
    for (int s1i = 0; s1i < static_cast<int>(map2.size()); ++s1i) {
      const Site& site1 = config.select_particle(particle_index1).site(s1i);
      const int site_type1 = site1.type();
      const std::vector<double>& map1 = smap[s0i][s1i];
      if (neighbor_criteria.is_accepted(map1[0], map1[1],
                                         site_type0, site_type1)) {
        if (frame) {
          frame->set_to_origin(dimen());
          for (int dim = 0; dim < dimen(); ++dim) {
            frame->set_coord(dim, -1.*map1[2 + dim]);
          }
        }
        return true;
      }
    }
  }
  return false;
}

void EnergyMapAll::check(const Configuration& config) const {
  if (!feasst::is_equal(map(), map_new(), NEAR_ZERO)) {
//    DEBUG(feasst_str(map()));
//    DEBUG(feasst_str(map_new()));
    ERROR("maps are not equal");
  }

  // check if ghost particles are in the neighbor list
  for (const int part : config.group_select(0).particle_indices()) {
    for (const std::shared_ptr<Select>& ghost : config.ghosts()) {
      for (int ghost_part : ghost->particle_indices()) {
        if (ghost_part < static_cast<int>(map()[part].size())) {
          for (int n_site = 0; n_site < static_cast<int>(map()[part][ghost_part].size()); ++n_site) {
            for (int g_site = 0; g_site < static_cast<int>(map()[part][ghost_part][n_site].size()); ++g_site) {
              if (map()[part][ghost_part][n_site][g_site][0] != 0) {
                INFO("existing particles: " << config.group_select(0).str());
                for (const std::shared_ptr<Select>& ghost2 : config.ghosts()) {
                  INFO("ghosts: " << ghost2->str());
                }
                FATAL("ghost particle " << ghost_part << " in map for p " << part << " s " << n_site);
              }
            }
          }
        }
      }
    }
  }
}

bool EnergyMapAll::is_cluster_changed(const NeighborCriteria& neighbor_criteria,
    const Select& select,
    const Configuration& config) const {
  for (int p1 : select.particle_indices()) {
    for (int p2 = 0; p2 < static_cast<int>(map()[p1].size()); ++p2) {
      if (is_cluster_(neighbor_criteria, map()[p1][p2], p1, p2, config) !=
          is_cluster_(neighbor_criteria, map_new()[p1][p2], p1, p2, config)) {
        return true;
      }
    }
  }
  return false;
}

typedef std::vector<std::vector<std::vector<std::vector<double> > > > vec4;

void EnergyMapAll::neighbors(
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
  DEBUG("site_type0 " << site_type0);
  DEBUG("target_particle " << target_particle);
  DEBUG("sz " << map().size());
  const vec4 * map4 = const_cast<vec4 * const>(&map()[target_particle]);
  if (new_map == 1) {
    map4 = const_cast<vec4 * const>(&map_new()[target_particle]);
  }
  for (int ipart = 0; ipart < static_cast<int>(map4->size()); ++ipart) {
    const Site& site1 = config.select_particle(ipart).site(given_site_index);
    const int site_type1 = site1.type();
    const std::vector<double> & map1 = (*map4)[ipart][target_site][given_site_index];
    DEBUG("site_type1 " << site_type1);
    if (neighbor_criteria.is_accepted(map1[0], map1[1], site_type0, site_type1)) {
      neighbors->add_site(ipart, given_site_index);
    }
  }
}

std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > * EnergyMapAll::map_() {
  return &((*data_.get_dble_6D())[0]);
}

std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > * EnergyMapAll::map_new_() {
  return &((*data_.get_dble_6D())[1]);
}

void EnergyMapAll::synchronize_(const EnergyMap& emap, const Select& perturbed) {
  for (int sel_index = 0; sel_index < perturbed.num_particles(); ++sel_index) {
    const int p1 = perturbed.particle_index(sel_index);
    for (const int s1 : perturbed.site_indices(sel_index)) {
      for (int p2 = 0; p2 < static_cast<int>(map()[p1].size()); ++p2) {
        for (int s2 = 0; s2 < static_cast<int>(map()[p1][p2][s1].size()); ++s2) {
//          ASSERT(emap.data().dble_6D().size() > 0, "err");
//          ASSERT(emap.data().dble_6D()[0].size() > p1, "err");
//          ASSERT(emap.data().dble_6D()[0][p1].size() > p2, "err");
//          ASSERT(emap.data().dble_6D()[0][p1][p2].size() > s1, "err");
//          ASSERT(emap.data().dble_6D()[0][p1][p2][s1].size() > s2, "err");
          (*map_())[p1][p2][s1][s2] = emap.data().dble_6D()[0][p1][p2][s1][s2];
          (*map_())[p2][p1][s2][s1] = emap.data().dble_6D()[0][p2][p1][s2][s1];
          (*map_new_())[p1][p2][s1][s2] = emap.data().dble_6D()[1][p1][p2][s1][s2];
          (*map_new_())[p2][p1][s2][s1] = emap.data().dble_6D()[1][p2][p1][s2][s1];
        }
      }
    }
  }
}

double EnergyMapAll::energy(const int part1_index, const int site1_index) const {
  double energy = 0.;
  for (const std::vector<std::vector<std::vector<double> > >& map3 : map()[part1_index]) {
    for (const std::vector<double>& map : map3[site1_index]) {
      energy += map[0];
    }
  }
  return energy;
}

}  // namespace feasst
