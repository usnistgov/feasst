#include "configuration/include/particle.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "utils/include/utils.h"
#include "math/include/utils_math.h"

namespace feasst {

void Particle::displace(const Position& displacement) {
  for (Site& site : sites_) {
    site.displace(displacement);
  }
}

void Particle::check() {
  for (Site site : sites_) {
    ASSERT(sites_[0].position().size() == site.position().size(), "size error");
  }
  properties().check();
}

void Particle::remove_non_unique_types() {
  std::vector<std::string> names = {"site", "bond", "angle", "dihedral"};
  for (std::string name : names) {
    int num = -1;
    if (name == "site") {
      num = num_sites();
    } else if (name == "bond") {
      num = num_bonds();
    } else if (name == "angle") {
      num = num_angles();
    } else if (name == "dihedral") {
      num = num_dihedrals();
    } else {
      ERROR("unrecognized");
    }
    std::vector<int> to_remove;
    std::vector<int> types_visited;
    for (int index = 0; index < num; ++index) {
      int type;
      if (name == "site") {
        type = sites_[index].type();
      } else if (name == "bond") {
        type = bonds_[index].type();
      } else if (name == "angle") {
        type = angles_[index].type();
      } else if (name == "dihedral") {
        type = dihedrals_[index].type();
      } else {
        FATAL("unrecognized");
      }
      if (!find_in_list(type, types_visited)) {
        types_visited.push_back(type);
      } else {
        to_remove.push_back(index);
      }
    }
    for (int index = num - 1; index >= 0; --index) {
      if (find_in_list(index, to_remove)) {
        if (name == "site") {
          sites_.erase(sites_.begin()  + index);
        } else if (name == "bond") {
          bonds_.erase(bonds_.begin()  + index);
        } else if (name == "angle") {
          angles_.erase(angles_.begin()  + index);
        } else if (name == "dihedral") {
          dihedrals_.erase(dihedrals_.begin()  + index);
        }
      }
    }
  }
}

void Particle::remove_site(const int index) {
  sites_.erase(sites_.begin() + index);
//  if (reference_sites_.size() != 0) {
//    reference_sites_.erase(sites_.begin() + index);
//  }
}

void Particle::increment_site_types(const int increment) {
  for (Site& site : sites_) {
    site.set_type(site.type() + increment);
  }
}

bool Particle::is_isotropic() {
  if (num_sites() <= 1) { return true; } else { return false; }
}

void Particle::replace_position(const Particle& particle) {
  ASSERT(num_sites() == particle.num_sites(), "size error");
  for (int site_index = 0; site_index < num_sites(); ++site_index) {
    sites_[site_index].set_position(particle.site(site_index).position());
  }
}

void Particle::replace_position(const int site_index,
                                const Position& replacement) {
  sites_[site_index].set_position(replacement);
}

// void Particle::add_or_set_site_property(const std::string name,
//                                         const double value) {
//   for (int site_index = 0; site_index < num_sites(); ++site_index) {
//     add_or_set_site_property(name, value, site_index);
//   }
// }

void Particle::add_bond_(const Bond& bond, const int index,
    std::vector<std::vector<int> > * list) {
  list->resize(num_sites());
  for (const int site : bond.site_indices()) {
    (*list)[site].push_back(index);
  }
}

void Particle::add_bond(const Bond& bond) {
  const int bond_index = static_cast<int>(bonds_.size());
  bonds_.push_back(bond);
  add_bond_(bond, bond_index, &bond_list_);

  bond_neighbors_.resize(num_sites());
  for (int site1_index = 0;
       site1_index < static_cast<int>(bond.site_indices().size()) - 1;
       ++site1_index) {
    const int site1 = bond.site_indices()[site1_index];
    for (int site2_index = site1_index + 1;
         site2_index < static_cast<int>(bond.site_indices().size());
         ++site2_index) {
      const int site2 = bond.site_indices()[site2_index];
      bond_neighbors_[site1].push_back(site2);
      bond_neighbors_[site2].push_back(site1);
    }
  }
}

void Particle::add_angle(const Angle& angle) {
  const int angle_index = static_cast<int>(angles_.size());
  angles_.push_back(angle);
  add_bond_(angle, angle_index, &angle_list_);

  // update neighbors
  angle_neighbors_.resize(num_sites());
  for (int site1_index = 0;
       site1_index < static_cast<int>(angle.site_indices().size()) - 2;
       ++site1_index) {
    const int site1 = angle.site_indices()[site1_index];
    for (int site2_index = site1_index + 1;
         site2_index < static_cast<int>(angle.site_indices().size()) - 1;
         ++site2_index) {
      const int site2 = angle.site_indices()[site2_index];
      for (int site3_index = site2_index + 1;
           site3_index < static_cast<int>(angle.site_indices().size());
           ++site3_index) {
        const int site3 = angle.site_indices()[site3_index];
        angle_neighbors_[site1].push_back({site2, site3});
        angle_neighbors_[site3].push_back({site2, site1});
      }
    }
  }
}

void Particle::add_dihedral(const Dihedral& dihedral) {
  const int dihedral_index = static_cast<int>(dihedrals_.size());
  dihedrals_.push_back(dihedral);
  add_bond_(dihedral, dihedral_index, &dihedral_list_);

  // update neighbors
  dihedral_neighbors_.resize(num_sites());
  for (int site1_index = 0;
       site1_index < static_cast<int>(dihedral.site_indices().size()) - 3;
       ++site1_index) {
    const int site1 = dihedral.site_indices()[site1_index];
    for (int site2_index = site1_index + 1;
         site2_index < static_cast<int>(dihedral.site_indices().size()) - 2;
         ++site2_index) {
      const int site2 = dihedral.site_indices()[site2_index];
      for (int site3_index = site2_index + 1;
           site3_index < static_cast<int>(dihedral.site_indices().size()) - 1;
           ++site3_index) {
        const int site3 = dihedral.site_indices()[site3_index];
        for (int site4_index = site3_index + 1;
             site4_index < static_cast<int>(dihedral.site_indices().size());
             ++site4_index) {
          const int site4 = dihedral.site_indices()[site4_index];
          dihedral_neighbors_[site1].push_back({site2, site3, site4});
          dihedral_neighbors_[site4].push_back({site3, site2, site1});
        }
      }
    }
  }
}

const Bond& Particle::bond(const int site_index1, const int site_index2) const {
  DEBUG("sites " << site_index1 << " " << site_index2);
  ASSERT(site_index1 < num_sites(),
    "site:" << site_index1 << " is > number of sites:" << num_sites());
  ASSERT(site_index2 < num_sites(),
    "site:" << site_index2 << " is > number of sites:" << num_sites());
  for (const int bond_index : bond_list_[site_index1]) {
    const Bond& bond = Particle::bond(bond_index);
    if ( (site_index1 == bond.site(0) && (site_index2 == bond.site(1))) ||
         (site_index1 == bond.site(1) && (site_index2 == bond.site(0)))) {
      return bond;
    }
  }
  FATAL("bond between " << site_index1 << " and " << site_index2 <<
        " not found.");
}

const Angle& Particle::angle(const int site_index1,
    const int site_index2,
    const int site_index3) const {
  DEBUG("sites " << site_index1 << " " << site_index2 << " " << site_index3);
  ASSERT(site_index1 < static_cast<int>(angle_list_.size()), "site_index1: " <<
    site_index1 << " >= number of angles: " << angle_list_.size() <<
    ". Perhaps an angle needs to be defined in this particle.");
  for (const int angle_index : angle_list_[site_index1]) {
    const Angle& angle = Particle::angle(angle_index);
    if (find_in_list(site_index1, angle.site_indices())) {
      if (find_in_list(site_index2, angle.site_indices())) {
        if (find_in_list(site_index3, angle.site_indices())) {
          return angle;
        }
      }
    }
  }
  FATAL("angle between " << site_index1 << "," << site_index2 << "," <<
    site_index3 << " not found.");
}

const Dihedral& Particle::dihedral(const int site_index1, const int site_index2,
    const int site_index3, const int site_index4) const {
  DEBUG("sites " << site_index1 << " " << site_index2 << " "
                 << site_index3 << " " << site_index4);
  for (const int dihedral_index : dihedral_list_[site_index1]) {
    const Dihedral& dihedral = Particle::dihedral(dihedral_index);
    if (find_in_list(site_index1, dihedral.site_indices())) {
      if (find_in_list(site_index2, dihedral.site_indices())) {
        if (find_in_list(site_index3, dihedral.site_indices())) {
          if (find_in_list(site_index4, dihedral.site_indices())) {
            return dihedral;
          }
        }
      }
    }
  }
  FATAL("dihedral between " << site_index1 << " - " << site_index2 << " - " <<
    site_index3 << " - " << site_index4 << " not found.");
}

void Particle::erase_bonds() {
  bonds_.clear();
  bond_list_.clear();
  bond_neighbors_.clear();
  angles_.clear();
  angle_list_.clear();
  dihedrals_.clear();
}

int Particle::num_sites_of_type(const int type) const {
  int num = 0;
  for (const Site& site : sites_) {
    if (site.type() == type) {
      ++num;
    }
  }
  return num;
}

void Particle::serialize(std::ostream& ostr) const {
  PropertiedEntity::serialize(ostr);
  feasst_serialize_version(366, ostr);
  feasst_serialize(type_, ostr);
  feasst_serialize_fstobj(sites_, ostr);
  feasst_serialize_fstobj(bonds_, ostr);
  feasst_serialize_fstobj(angles_, ostr);
  feasst_serialize_fstobj(dihedrals_, ostr);
  feasst_serialize(bond_list_, ostr);
  feasst_serialize(bond_neighbors_, ostr);
  feasst_serialize(angle_list_, ostr);
  feasst_serialize(angle_neighbors_, ostr);
  feasst_serialize(dihedral_list_, ostr);
  feasst_serialize(dihedral_neighbors_, ostr);
}

Particle::Particle(std::istream& istr) : PropertiedEntity(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 365 && version <= 366, "version mismatch: " << version);
  feasst_deserialize(&type_, istr);
  feasst_deserialize_fstobj(&sites_, istr);
  feasst_deserialize_fstobj(&bonds_, istr);
  feasst_deserialize_fstobj(&angles_, istr);
  if (version >= 336) {
    feasst_deserialize_fstobj(&dihedrals_, istr);
  }
  feasst_deserialize(&bond_list_, istr);
  feasst_deserialize(&bond_neighbors_, istr);
  feasst_deserialize(&angle_list_, istr);
  if (version >= 336) {
    feasst_deserialize(&angle_neighbors_, istr);
    feasst_deserialize(&dihedral_list_, istr);
    feasst_deserialize(&dihedral_neighbors_, istr);
  }
}

const std::vector<int>& Particle::bond_neighbors(const int site) const {
  ASSERT(site < num_sites(), "site:" << site << " > num_sites: " <<
         num_sites());
  if (site >= static_cast<int>(bond_neighbors_.size())) {
    return empty_;
  }
  return bond_neighbors_[site];
}

const std::vector<std::vector<int> >& Particle::angle_neighbors(
    const int site) const {
  ASSERT(site < num_sites(), "site:" << site << " > num_sites: " <<
         num_sites());
  if (site >= static_cast<int>(angle_neighbors_.size())) {
    return empty2d_;
  }
  return angle_neighbors_[site];
}

const std::vector<std::vector<int> >& Particle::dihedral_neighbors(
    const int site) const {
  ASSERT(site < num_sites(), "site:" << site << " > num_sites: " <<
         num_sites());
  if (site >= static_cast<int>(dihedral_neighbors_.size())) {
    return empty2d_;
  }
  return dihedral_neighbors_[site];
}

double Particle::max_distance() const {
  double max_dist = 0.;
  for (const Site& site : sites_) {
    const double dist_sq = site.position().squared_distance();
    if (dist_sq > max_dist) {
      max_dist = dist_sq;
    }
  }
  return std::sqrt(max_dist);
}

const Site& Particle::site(const int index) const {
  ASSERT(index < num_sites(),
    "index: " << index << " >= num_sites: " << num_sites());
  return sites_[index];
}

void Particle::set_site(const int index, const Site& site) {
  ASSERT(index < num_sites(),
    "index: " << index << " >= num_sites: " << num_sites());
  sites_[index] = site;
}

void Particle::clear_names() {
  for (Site& site : sites_) {
    site.set_name("");
  }
  for (Bond& bond : bonds_) {
    bond.set_name("");
  }
  for (Angle& angle : angles_) {
    angle.set_name("");
  }
  for (Dihedral& dihedral : dihedrals_) {
    dihedral.set_name("");
  }
}

}  // namespace feasst
