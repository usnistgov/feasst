#include "cluster/include/energy_map_all.h"

namespace feasst {

class MapEnergyMapAll {
 public:
  MapEnergyMapAll() {
    EnergyMapAll().deserialize_map()["EnergyMapAll"] =
      MakeEnergyMapAll();
  }
};

static MapEnergyMapAll mapper_ = MapEnergyMapAll();

EnergyMapAll::EnergyMapAll(std::istream& istr) : EnergyMap(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 210, "mismatch:" << version);
  feasst_deserialize(&map_, istr);
}

void EnergyMapAll::serialize_energy_map_all_(std::ostream& ostr) const {
  feasst_serialize_version(210, ostr);
  feasst_serialize(map_, ostr);
}

void EnergyMapAll::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_energy_map_all_(ostr);
}

void EnergyMapAll::resize_(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index) {
  const int part_max_old = static_cast<int>(map_.size());
  if (part1_index >= part_max_old || part2_index >= part_max_old) {
    const int part_max = std::max(part1_index, part2_index) + 1;
    map_.resize(part_max);
    map_new_.resize(part_max);
    for (int i = 0; i < part_max; ++i) {
      map_[i].resize(part_max);
      map_new_[i].resize(part_max);
      for (int j = 0; j < part_max; ++j) {
        ASSERT(site_max() != 0, "wasn't precomputed");
        map_[i][j].resize(site_max());
        map_new_[i][j].resize(site_max());
        for (int k = 0; k < site_max(); ++k) {
          map_[i][j][k].resize(site_max());
          map_new_[i][j][k].resize(site_max());
          if (i >= part_max_old || j >= part_max_old) {
            for (int l = 0; l < site_max(); ++l) {
              ASSERT(dimen() != -1, "wasnt precomputed");
              map_[i][j][k][l] = std::vector<double>(2 + dimen(), default_value());
              map_new_[i][j][k][l] = map_[i][j][k][l];
            }
          }
        }
      }
    }
  }
}

void EnergyMapAll::revert(const Select& select) {
  const int pmax = select.particle_indices().back();
  const int smax = select.site_indices().back().back();
  resize_(pmax, smax, pmax, smax);
  for (int sel_index = 0; sel_index < select.num_particles(); ++sel_index) {
    const int p1 = select.particle_index(sel_index);
    for (const int s1 : select.site_indices(sel_index)) {
      for (int p2 = 0; p2 < static_cast<int>(map_[p1].size()); ++p2) {
        for (int s2 = 0; s2 < static_cast<int>(map_[p1][p2][s1].size()); ++s2) {
          if (select.trial_state() == 3) {
            // revert addition
            map_[p1][p2][s1][s2] = std::vector<double>(2 + dimen(), default_value());
            map_[p2][p1][s2][s1] = std::vector<double>(2 + dimen(), default_value());
            map_new_[p1][p2][s1][s2] = std::vector<double>(2 + dimen(), default_value());
            map_new_[p2][p1][s2][s1] = std::vector<double>(2 + dimen(), default_value());
          } else {
            map_new_[p1][p2][s1][s2] = map_[p1][p2][s1][s2];
            map_new_[p2][p1][s2][s1] = map_[p2][p1][s2][s1];
          }
        }
      }
    }
  }
}

void EnergyMapAll::finalize(const Select& select) {
  if (select.num_particles() > 0) {
    const int pmax = select.particle_indices().back();
    const int smax = select.site_indices().back().back();
    resize_(pmax, smax, pmax, smax);
    for (int sel_index = 0; sel_index < select.num_particles(); ++sel_index) {
      const int p1 = select.particle_index(sel_index);
      for (const int s1 : select.site_indices(sel_index)) {
        for (int p2 = 0; p2 < static_cast<int>(map_[p1].size()); ++p2) {
          for (int s2 = 0; s2 < static_cast<int>(map_[p1][p2][s1].size()); ++s2) {
            TRACE("p1s1p2s2 " << p1 << " " << s1 << " " << p2 << " " << s2);
            if (select.trial_state() == 2) {
              // finalize removal
              map_[p1][p2][s1][s2] = std::vector<double>(2 + dimen(), default_value());
              map_[p2][p1][s2][s1] = std::vector<double>(2 + dimen(), default_value());
              map_new_[p1][p2][s1][s2] = std::vector<double>(2 + dimen(), default_value());
              map_new_[p2][p1][s2][s1] = std::vector<double>(2 + dimen(), default_value());
            } else {
              map_[p1][p2][s1][s2] = map_new_[p1][p2][s1][s2];
              map_[p2][p1][s2][s1] = map_new_[p2][p1][s2][s1];
            }
          }
        }
      }
    }
  }
}

void EnergyMapAll::select_cluster(const ClusterCriteria * cluster_criteria,
                                  const Configuration& config,
                                  const int particle_node,
                                  SelectPosition * cluster,
                                  const Position& frame_of_reference) const {
  DEBUG("particle_node " << particle_node);
  DEBUG("map size " << map_.size());
  for (int part2_index = 0;
       part2_index < static_cast<int>(map_[particle_node].size());
       ++part2_index) {
    DEBUG("part2_index " << part2_index);
    // if part2 isn't already in the cluster
    // and part2 satistifies cluster criteria,
    // then recurively add part2 as a new node
    if (!find_in_list(part2_index, cluster->particle_indices())) {
      Position frame;
      if (is_cluster_(cluster_criteria, &map_[particle_node][part2_index], &frame)) {
        DEBUG("FOR " << frame_of_reference.str());
        frame.add(frame_of_reference);
        //frame.multiply(-1.);
        const Particle& part = config.select_particle(part2_index);
        cluster->add_particle(part, part2_index);
        cluster->load_positions_of_last(part, frame);
        DEBUG("frame: " << frame.str());
        DEBUG("added: " << cluster->site_positions().back().back().str());
        select_cluster(cluster_criteria, config, part2_index, cluster, frame);
      }
    }
  }
}

bool EnergyMapAll::is_cluster_(
    const ClusterCriteria * cluster_criteria,
    const std::vector<std::vector<std::vector<double> > > * smap,
    Position * frame) const {
  for (const std::vector<std::vector<double> > & map2 : *smap) {
    for (const std::vector<double>& map1 : map2) {
      if (cluster_criteria->is_accepted(map1[0], map1[1])) {
        if (frame != NULL) {
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

void EnergyMapAll::check() const {
  if (!is_equal(map_, map_new_, NEAR_ZERO)) {
    INFO(feasst_str(map_));
    INFO(feasst_str(map_new_));
    ERROR("maps are not equal");
  }
}

bool EnergyMapAll::is_cluster_changed(const ClusterCriteria * cluster_criteria,
    const Select& select) const {
  for (int p1 : select.particle_indices()) {
    for (int p2 = 0; p2 < static_cast<int>(map_[p1].size()); ++p2) {
      if (is_cluster_(cluster_criteria, &map_[p1][p2]) !=
          is_cluster_(cluster_criteria, &map_new_[p1][p2])) {
        return true;
      }
    }
  }
  return false;
}

}  // namespace feasst
