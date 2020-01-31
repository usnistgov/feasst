#include "system/include/energy_map_all.h"

namespace feasst {

class MapEnergyMapAll {
 public:
  MapEnergyMapAll() {
    EnergyMapAll().deserialize_map()["EnergyMapAll"] =
      MakeEnergyMapAll();
  }
};

static MapEnergyMapAll mapper_ = MapEnergyMapAll();

EnergyMapAll::EnergyMapAll(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 210, "mismatch:" << version);
  feasst_deserialize(&map_, istr);
  feasst_deserialize(&site_max_, istr);
  feasst_deserialize(&dimen_, istr);
}

void EnergyMapAll::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(210, ostr);
  feasst_serialize(map_, ostr);
  feasst_serialize(site_max_, ostr);
  feasst_serialize(dimen_, ostr);
}

void EnergyMapAll::precompute(Configuration * config) {
  site_max_ = config->max_sites_in_any_particle();
  DEBUG("site_max_ " << site_max_);
  dimen_ = config->dimension();
}

void EnergyMapAll::clear(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index) {
  resize_(part1_index, site1_index, part2_index, site2_index);
  map_[part1_index][part2_index][site1_index][site2_index] =
    std::vector<double>(2 + dimen_, default_value());
  map_[part2_index][part1_index][site2_index][site1_index] =
    std::vector<double>(2 + dimen_, default_value());
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
    map_old_.resize(part_max);
    for (int i = 0; i < part_max; ++i) {
      map_[i].resize(part_max);
      map_old_[i].resize(part_max);
      for (int j = 0; j < part_max; ++j) {
        ASSERT(site_max_ != 0, "wasn't precomputed");
        map_[i][j].resize(site_max_);
        map_old_[i][j].resize(site_max_);
        for (int k = 0; k < site_max_; ++k) {
          map_[i][j][k].resize(site_max_);
          map_old_[i][j][k].resize(site_max_);
          if (i >= part_max_old || j >= part_max_old) {
            for (int l = 0; l < site_max_; ++l) {
              ASSERT(dimen_ != -1, "wasnt precomputed");
              map_[i][j][k][l] = std::vector<double>(2 + dimen_, default_value());
              //clear(i, j, k, l);
              map_old_[i][j][k][l] = map_[i][j][k][l];
            }
          }
        }
      }
    }
  }
}

double EnergyMapAll::update(
    const double energy,
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const double squared_distance,
    const Position * pbc) {
  resize_(part1_index, site1_index, part2_index, site2_index);
  DEBUG("map: energy: " << energy << " p1p2s1s2: " << part1_index << " " << part2_index << " " << site1_index << " " << site2_index);
  DEBUG(map_.size() << " " << map_[0].size() << " " << map_[0][0].size() << " " << map_[0][0][0].size() << " " << map_[0][0][0][0].size());
  std::vector<double> * smap1 =
    &map_[part1_index][part2_index][site1_index][site2_index];
  (*smap1)[0] = energy;
  (*smap1)[1] = squared_distance;
  std::vector<double> * smap2 =
    &map_[part2_index][part1_index][site2_index][site1_index];
  (*smap2)[0] = energy;
  (*smap2)[1] = squared_distance;
  if (pbc->dimension() > 0) {
    for (int dim = 0; dim < dimen_; ++dim) {
      DEBUG("dim " << dim << " sz " << map_[part1_index][part2_index][site1_index][site2_index].size());
      DEBUG("pbc sz " << pbc->size());
      (*smap1)[2 + dim] = pbc->coord(dim);
      (*smap2)[2 + dim] = -1.*pbc->coord(dim);
    }
  }
  return energy;
}

double EnergyMapAll::query(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index) {
  return map_[part1_index][part2_index][site1_index][site2_index][0];
}
//void EnergyMapAll::revert() {
//  DEBUG("hi");
//}

void EnergyMapAll::remove_particles(const Select& select) {
  DEBUG("sel: " << select.str());
  const int pmax = select.particle_indices().back();
  const int smax = select.site_indices().back().back();
  resize_(pmax, smax, pmax, smax);
  DEBUG("map size: " << map_.size());
  DEBUG("map size: " << map_[0].size());
  DEBUG("map size: " << map_[0][0].size());
  DEBUG("map size: " << map_[0][0][0].size());
  DEBUG("map: " << feasst_str(map_));
  for (int sel_index = 0; sel_index < select.num_particles(); ++sel_index) {
    const int pindex = select.particle_index(sel_index);
    for (const int sindex : select.site_indices(sel_index)) {
      for (int p1 = 0; p1 < static_cast<int>(map_[pindex].size()); ++p1) {
        for (int s1 = 0; s1 < static_cast<int>(map_[pindex][p1][sindex].size()); ++s1) {
          map_[pindex][p1][sindex][s1] = std::vector<double>(2 + dimen_, default_value());
          map_[p1][pindex][s1][sindex] = std::vector<double>(2 + dimen_, default_value());
        }
      }
    }
  }
}

void EnergyMapAll::prep_for_revert(const Select& select) {
  const int pmax = select.particle_indices().back();
  const int smax = select.site_indices().back().back();
  resize_(pmax, smax, pmax, smax);
  for (int sel_index = 0; sel_index < select.num_particles(); ++sel_index) {
    const int pindex = select.particle_index(sel_index);
    for (const int sindex : select.site_indices(sel_index)) {
      for (int p1 = 0; p1 < static_cast<int>(map_[pindex].size()); ++p1) {
        for (int s1 = 0; s1 < static_cast<int>(map_[pindex][p1][sindex].size()); ++s1) {
          map_old_[pindex][p1][sindex][s1] = map_[pindex][p1][sindex][s1];
          map_old_[p1][pindex][s1][sindex] = map_[p1][pindex][s1][sindex];
        }
      }
    }
  }
}

void EnergyMapAll::revert(const Select& select) {
  for (int sel_index = 0; sel_index < select.num_particles(); ++sel_index) {
    const int pindex = select.particle_index(sel_index);
    for (const int sindex : select.site_indices(sel_index)) {
      for (int p1 = 0; p1 < static_cast<int>(map_[pindex].size()); ++p1) {
        for (int s1 = 0; s1 < static_cast<int>(map_[pindex][p1][sindex].size()); ++s1) {
          map_[pindex][p1][sindex][s1] = map_old_[pindex][p1][sindex][s1];
          map_[p1][pindex][s1][sindex] = map_old_[p1][pindex][s1][sindex];
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
    const std::vector<std::vector<std::vector<double> > >& smap =
      map_[particle_node][part2_index];
    if (!find_in_list(part2_index, cluster->particle_indices())) {
      Position frame;
      if (is_cluster_(cluster_criteria, smap, &frame)) {
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
    const std::vector<std::vector<std::vector<double> > >& smap,
    Position * frame) const {
  for (const std::vector<std::vector<double> > & map2 : smap) {
    for (const std::vector<double>& map1 : map2) {
      if (cluster_criteria->is_accepted(map1)) {
        frame->set_to_origin(dimen_);
        for (int dim = 0; dim < dimen_; ++dim) {
          frame->set_coord(dim, -1.*map1[2 + dim]);
        }
        return true;
      }
    }
  }
  return false;
}

double EnergyMapAll::total_energy() const {
  double en = 0;
  //DEBUG(map_.size() << " " << map_[0].size() << " " << map_[0][0].size() << " " << map_[0][0][0].size() << " " << map_[0][0][0][0].size());
  for (const std::vector<std::vector<std::vector<std::vector<double> > > >& vec4 : map_) {
    for (const std::vector<std::vector<std::vector<double> > >& vec3 : vec4) {
      for (const std::vector<std::vector<double> >& vec2 : vec3) {
        for (const std::vector<double>& vec1 : vec2) {
          DEBUG(vec1.size());
          en += vec1[0];
        }
      }
    }
  }
  return 0.5*en;
}

}  // namespace feasst
