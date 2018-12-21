
#include "core/include/group.h"
#include "core/include/utils.h"

namespace feasst {

Group::Group() { set_dynamic(); }

bool Group::empty() const {
  if ( (particle_types_.size() == 0) &&
       (site_types_.size() == 0) ) {
    return true;
  }
  return false;
}

bool Group::is_in(const Site& site) const {
  if (site_types_.size() == 0) {
    return true;
  }
  if (find_in_list(site.type(), site_types_)) {
    return true;
  }
  return false;
}

bool Group::is_in(const Particle& particle) const {
  if (particle_types_.size() == 0) {
    return true;
  }
  if (find_in_list(particle.type(), particle_types_)) {
    return true;
  }
  return false;
}

void Group::remove_sites(Particle * particle) const {
//                         std::vector<int> * full_to_partial,
//                         std::vector<int> * partial_to_full) const {
//  // compute site mappings (if requested) before removal
//  if (full_to_partial != NULL && partial_to_full != NULL) {
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

///// HWH depreciate after partials are removed
//Particle Group::remove_sites(const Particle& particle) const {
////                             std::vector<int> * full_to_partial,
////                             std::vector<int> * partial_to_full) const {
//  Particle filtered(particle);
//  remove_sites(&filtered);//, full_to_partial, partial_to_full);
//  return filtered;
//}

std::vector<int> Group::site_indices(const Particle& particle) const {
  std::vector<int> indices;
  for (int index = 0; index < particle.num_sites(); ++index) {
    if (is_in(particle.site(index))) {
      indices.push_back(index);
    }
  }
  return indices;
}

}  // namespace feasst
