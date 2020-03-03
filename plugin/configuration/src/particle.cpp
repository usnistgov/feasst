#include "configuration/include/particle.h"
#include "utils/include/debug.h"
#include "utils/include/utils_io.h"
#include "utils/include/utils.h"
#include "math/include/utils_math.h"

namespace feasst {

void Particle::displace(const Position& displacement) {
  add_position(displacement);
  for (Site& site : sites_) {
    site.displace(displacement);
  }
}

void Particle::check() {
  for (Site site : sites_) {
    ASSERT(position().size() == site.position().size(), "size error");
  }
  properties().check();
}

void Particle::remove_non_unique_types() {
  std::vector<std::string> names = {"site", "bond"};
  for (std::string name : names) {
    int num = -1;
    if (name == "site") {
      num = num_sites();
    } else if (name == "bond") {
      num = num_bonds();
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
      } else {
        ERROR("unrecognized");
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
  set_position(particle.position());
  ASSERT(num_sites() == particle.num_sites(), "size error");
  for (int site_index = 0; site_index < num_sites(); ++site_index) {
    sites_[site_index].set_position(particle.site(site_index).position());
  }
}

void Particle::replace_position(const int site_index,
                                const Position& replacement) {
  sites_[site_index].set_position(replacement);
}

//void Particle::add_or_set_site_property(const std::string name,
//                                        const double value) {
//  for (int site_index = 0; site_index < num_sites(); ++site_index) {
//    TRACE("name " << name << " value " << value << " site index " << site_index);
//    add_or_set_site_property(name, value, site_index);
//  }
//}

Position Particle::average_site_position() const {
  Position center = site(0).position();
  center.multiply(0);
  for (const Site& site : sites()) {
    center.add(site.position());
  }
  center.divide(num_sites());
  return center;
}

void Particle::set_position_as_center() {
  set_position(average_site_position());
}

void Particle::resize_list_(std::vector<std::vector<int> > * list) {
  if (static_cast<int>(list->size()) < num_sites()) {
    list->resize(num_sites());
  }
}

void Particle::add_bond_(const Bond& bond, const int index,
    std::vector<std::vector<int> > * list) {
  resize_list_(list);
  for (const int site : bond.site_indices()) {
    (*list)[site].push_back(index);
  }
}

void Particle::add_bond_neighbor_(const Bond& bond,
    std::vector<std::vector<int> > * list) {
  resize_list_(list);
  for (int site1_index = 0;
       site1_index < static_cast<int>(bond.site_indices().size()) - 1;
       ++site1_index) {
    const int site1 = bond.site_indices()[site1_index];
    for (int site2_index = site1_index + 1;
         site2_index < static_cast<int>(bond.site_indices().size());
         ++site2_index) {
      const int site2 = bond.site_indices()[site2_index];
      (*list)[site1].push_back(site2);
      (*list)[site2].push_back(site1);
    }
  }
}

void Particle::add_bond(const Bond& bond) {
  const int bond_index = static_cast<int>(bonds_.size());
  bonds_.push_back(bond);
  add_bond_(bond, bond_index, &bond_list_);
  add_bond_neighbor_(bond, &bond_neighbor_);
}

void Particle::add_angle(const Angle& angle) {
  const int angle_index = static_cast<int>(angles_.size());
  angles_.push_back(angle);
  add_bond_(angle, angle_index, &angle_list_);
}

const Bond& Particle::bond(const int site_index1, const int site_index2) const {
  DEBUG("sites " << site_index1 << " " << site_index2);
  for (const int bond_index : bond_list_[site_index1]) {
    const Bond& bond = Particle::bond(bond_index);
    if ( (site_index1 == bond.site(0) and (site_index2 == bond.site(1))) or
         (site_index1 == bond.site(1) and (site_index2 == bond.site(0)))) {
      return bond;
    }
  }
  FATAL("bond between " << site_index1 << " and " << site_index2 << " not found.");
}

const Angle& Particle::angle(const int site_index1,
    const int site_index2,
    const int site_index3) const {
  DEBUG("sites " << site_index1 << " " << site_index2 << " " << site_index3);
  for (const int angle_index : angle_list_[site_index1]) {
    const Angle& angle = Particle::angle(angle_index);
    if ( ( (site_index1 == angle.site(0)) and
           (site_index2 == angle.site(1)) and
           (site_index3 == angle.site(2)) ) or
         ( (site_index1 == angle.site(2)) and
           (site_index2 == angle.site(1)) and
           (site_index3 == angle.site(0)) ) ) {
      return angle;
    }
  }
  FATAL("angle between " << site_index1 << " - " << site_index2 << " - " <<
    site_index3 << " not found.");
}

void Particle::erase_bonds() {
  bonds_.clear();
  bond_list_.clear();
  bond_neighbor_.clear();
  angles_.clear();
  angle_list_.clear();
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
  TypedEntity::serialize(ostr);
  SpatialEntity::serialize(ostr);
  feasst_serialize_version(365, ostr);
  feasst_serialize_fstobj(sites_, ostr);
  feasst_serialize_fstobj(bonds_, ostr);
  feasst_serialize_fstobj(angles_, ostr);
  feasst_serialize(bond_list_, ostr);
  feasst_serialize(bond_neighbor_, ostr);
  feasst_serialize(angle_list_, ostr);
}

Particle::Particle(std::istream& istr)
  : PropertiedEntity(istr),
    TypedEntity(istr),
    SpatialEntity(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 365, "version mismatch: " << version);
  feasst_deserialize_fstobj(&sites_, istr);
  feasst_deserialize_fstobj(&bonds_, istr);
  feasst_deserialize_fstobj(&angles_, istr);
  feasst_deserialize(&bond_list_, istr);
  feasst_deserialize(&bond_neighbor_, istr);
  feasst_deserialize(&angle_list_, istr);
}

}  // namespace feasst
