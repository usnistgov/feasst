
#include "utils/include/arguments.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "configuration/include/particle.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/group.h"

namespace feasst {

Group::Group(argtype * args) : PropertiedEntity() {
  std::string pre = str("prepend", args, "");
  if (!pre.empty()) {
    pre += "_";
  }

  std::string start;
  // if only one site type, allow the subscript to be dropped
  start.assign(pre+"site_type");
  if (used(start, *args)) {
    for (const std::string& st : split(feasst::str(start, args), ',')) {
      site_type_names_.push_back(st);
      try {
        site_types_.push_back(std::stoi(site_type_names_.back()));
      } catch (...) {
        site_types_.push_back(-1);
      }
    }
  } else {
    int type = site_type_names_.size();
    std::stringstream key;
    key << start << type;
    while (used(key.str(), *args)) {
      WARN("site_type[i] is deprecated. Use comma-separated input.");
      site_type_names_.push_back(str(key.str(), args));
      try {
        site_types_.push_back(std::stoi(site_type_names_.back()));
      } catch (...) {
        site_types_.push_back(-1);
      }
      ++type;
      ASSERT(type < 1e8, "type(" << type << ") is very high. Infinite loop?");
      key.str("");
      key << start << type;
    }
  }

  start.assign(pre+"particle_type");
  if (used(start, *args)) {
    for (const std::string& st : split(feasst::str(start, args), ',')) {
      particle_type_names_.push_back(st);
      try {
        particle_types_.push_back(std::stoi(particle_type_names_.back()));
      } catch (...) {
        particle_types_.push_back(-1);
      }
    }
  } else {
    int type = particle_type_names_.size();
    std::stringstream key;
    key << start << type;
    while (used(key.str(), *args)) {
      WARN("particle_type[i] is deprecated. Use comma-separated input.");
      particle_type_names_.push_back(str(key.str(), args));
      try {
        particle_types_.push_back(std::stoi(particle_type_names_.back()));
      } catch (...) {
        particle_types_.push_back(-1);
      }
      ++type;
      ASSERT(type < 1e8, "type(" << type << ") is very high. Infinite loop?");
      key.str("");
      key << start << type;
    }
  }

  // if only one particle index, allow drop the subscript
  start.assign(pre+"particle_index");
  if (used(start, *args)) {
    for (const std::string& st : split(feasst::str(start, args), ',')) {
      particle_indices_.push_back(str_to_int(st));
    }
  } else {
    int index = particle_indices_.size();
    std::stringstream key;
    key << start << index;
    while (used(key.str(), *args)) {
      WARN("particle_type[i] is deprecated. Use comma-separated input.");
      particle_indices_.push_back(integer(key.str(), args));
      ++index;
      ASSERT(index < 1e8,
             "index(" << index << ") is very high. Infinite loop?");
      key.str("");
      key << start << index;
    }
  }

  dynamic_ = boolean(pre+"dynamic", args, true);
  spatial_ = boolean(pre+"spatial", args, false);
  ASSERT(!spatial_, "spatial groups are not implemented");
}
Group::Group(argtype args) : Group(&args) {
  feasst_check_all_used(args);
}

void Group::name_to_index(const ParticleFactory& unique_types) {
  site_types_.clear();
  for (const std::string& stn : site_type_names_) {
    site_types_.push_back(unique_types.site_type_name_to_index(stn));
  }
  particle_types_.clear();
  for (const std::string& ptn : particle_type_names_) {
    particle_types_.push_back(unique_types.name_to_index(ptn));
  }
}

bool Group::is_empty() const {
  if ( (particle_types_.size() == 0) && (site_types_.size() == 0) ) {
    return true;
  }
  return false;
}

bool Group::is_in(const Site& site) const {
  bool type = false;
  if ((site_types_.size() == 0) ||
      find_in_list(site.type(), site_types_)) {
    type = true;
  }
  return type;
}

bool Group::is_in(const Particle& particle, const int particle_index) const {
  bool type = false;
  if ((particle_types_.size() == 0) ||
      find_in_list(particle.type(), particle_types_)) {
    type = true;
  }
  bool index = false;
  if ((particle_indices_.size() == 0) ||
      find_in_list(particle_index, particle_indices_)) {
    index = true;
  }
  return (type && index);
}

void Group::remove_sites(Particle * particle) const {
//                         std::vector<int> * full_to_partial,
//                         std::vector<int> * partial_to_full) const {
//  // compute site mappings (if requested) before removal
//  if (full_to_partial != NULL and partial_to_full != NULL) {
//    int partial_site = 0;
//    for (int index = 0; index < particle->num_sites(); ++index) {
//      if (is_in(particle->site(index))) {
//        (*full_to_partial).push_back(partial_site);
//        ++partial_site;
//        (*partial_to_full).push_back(index);
//      } else {
//        (*full_to_partial).push_back(-1);
//      }
//    }
//  }

  // loop backward to remove sites
  for (int index = particle->num_sites() - 1; index >= 0; --index) {
    if (!is_in(particle->site(index))) {
      particle->remove_site(index);
    }
  }
}

// Particle Group::remove_sites(const Particle& particle) const {
// //                             std::vector<int> * full_to_partial,
// //                             std::vector<int> * partial_to_full) const {
//   Particle filtered(particle);
//   remove_sites(&filtered);//, full_to_partial, partial_to_full);
//   return filtered;
// }

std::vector<int> Group::site_indices(const Particle& particle) const {
  std::vector<int> indices;
  for (int index = 0; index < particle.num_sites(); ++index) {
    if (is_in(particle.site(index))) {
      indices.push_back(index);
    }
  }
  return indices;
}

void Group::serialize(std::ostream& ostr) const {
  PropertiedEntity::serialize(ostr);
  feasst_serialize_version(1036, ostr);
  feasst_serialize(site_types_, ostr);
  feasst_serialize(site_type_names_, ostr);
  feasst_serialize(particle_type_names_, ostr);
  feasst_serialize(particle_types_, ostr);
  // feasst_serialize(site_indices_, ostr);
  feasst_serialize(particle_indices_, ostr);
  ostr << dynamic_ << " " << spatial_ << " ";
}

Group::Group(std::istream& istr) : PropertiedEntity(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 1035 && version <= 1036, "unrecognized version: " << version);
  feasst_deserialize(&site_types_, istr);
  if (version >= 1036) {
    feasst_deserialize(&site_type_names_, istr);
    feasst_deserialize(&particle_type_names_, istr);
  }
  feasst_deserialize(&particle_types_, istr);
  // feasst_deserialize(&site_indices_, istr);
  feasst_deserialize(&particle_indices_, istr);
  istr >> dynamic_ >> spatial_;
}

}  // namespace feasst
