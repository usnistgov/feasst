
#include <sstream>
#include <algorithm>
#include "core/include/select.h"
#include "core/include/utils.h"
#include "core/include/utils_math.h"
#include "core/include/utils_io.h"
#include "core/include/debug.h"

namespace feasst {

//void Select::add_particle(const Particles& particles,
//                             const int particle_index) {
//  const Particle& part = particles.particle(particle_index);
//  add_particle(part, particle_index);
//}

void Select::add(const Select& select) {
  TRACE("adding " << select.str() << " to " << this->str());
  for (int select_index = 0;
       select_index < select.num_particles();
       ++select_index) {
    const int particle_index = select.particle_index(select_index);
    const std::vector<int>& site_indices = select.site_indices(select_index);
    add_sites(particle_index, site_indices);
  }
}

void Select::remove(const Select& select) {
  for (int select_index = 0; select_index < select.num_particles(); ++select_index) {
    remove_sites(select.particle_index(select_index),
                 select.site_indices(select_index));
  }
}

void Select::add_particle(const Particle& particle,
    const int particle_index) {
  std::vector<int> sites;
  for (int index = 0; index < particle.num_sites(); ++index) {
    sites.push_back(index);
  }
  particle_indices_.push_back(particle_index);
  site_indices_.push_back(sites);
}

std::string Select::str() const {
  std::stringstream ss;
  for (int index = 0; index < num_particles(); ++index) {
    ss << particle_indices_[index] << ":{"
       << feasst_str(site_indices_[index]) << "}, ";
  }
  return ss.str();
}

bool Select::is_empty() const {
  if (num_particles() == 0) {
    return true;
  }
  return false;
}

//int Select::random_particle_index(const Particles& particles) {
//  ASSERT(particles.num() > 0, "size error");
//  return static_cast<int>(particles.num() * random_.uniform());
//}

int Select::num_sites() const {
  return num_elements(site_indices_);
}

Select Select::random_particle() {
  Select select;
  if (num_particles() > 0) {
    int select_index;
    const int particle_index = random_.element(particle_indices_,
                                               &select_index);
    DEBUG("num " << num_particles());
    DEBUG("index " << particle_index);
    const std::vector<int> sites = site_indices_[select_index];
    select.add_particle(particle_index, sites);
  }
  return select;
}

/// HWH depreciate
void Select::reduce_to_random_particle() {
  if (num_particles() <=0) {
    return;
  }
  const int particle_index = random_.element(particle_indices_);
  const std::vector<int> sites = site_indices_[particle_index];
  clear();
  particle_indices_.push_back(particle_index);
  site_indices_.push_back(sites);
}

void Select::add_site(const int particle_index, const int site_index) {
  int index;
  if (find_in_list(particle_index, particle_indices(), &index)) {
    site_indices_[index].push_back(site_index);
  } else {
    particle_indices_.push_back(particle_index);
    site_indices_.push_back({site_index});
  }
}

void Select::add_sites(const int particle_index,
                       const std::vector<int> site_indices) {
  int index;
  if (find_in_list(particle_index, particle_indices(), &index)) {
    TRACE(particle_index << " | " << feasst_str(site_indices_[index]) << " + " << feasst_str(site_indices));
    site_indices_[index] = fst_union(site_indices_[index], site_indices);
  } else {
    particle_indices_.push_back(particle_index);
    site_indices_.push_back(site_indices);
  }
}

void Select::remove_sites(const int particle_index,
                          const std::vector<int> site_indices) {
  int index;
  if (find_in_list(particle_index, particle_indices(), &index)) {
    site_indices_[index] = fst_difference(site_indices_[index], site_indices);
    // completely remove particle if no more sites remain.
    if (site_indices_[index].size() == 0) {
      remove_particle_(index);
    }
  }
}

//void Select::add_last_particle(const Particles& particles) {
//  add_particle(particles, particles.num() - 1);
//}

void Select::remove_last_particle() {
  particle_indices_.pop_back();
  site_indices_.pop_back();
}

void Select::remove_particle(const int particle_index) {
  int index;
  if (find_in_list(particle_index, particle_indices_, &index)) {
    remove_particle_(index);
  }
}

void Select::check_size() const {
  ASSERT(particle_indices_.size() == site_indices_.size(), "size error");
}

void Select::add_particle(const int particle_index, std::vector<int> site_indices) {
  if (site_indices.size() > 0) {
    particle_indices_.push_back(particle_index);
    site_indices_.push_back(site_indices);
  }
}

bool Select::is_equivalent(const Select& select) {
  bool equal = true;
  if (particle_indices_ != select.particle_indices()) {
    equal = false;
  }
  if (site_indices_ != select.site_indices()) {
    equal = false;
  }
  if (unique_id_ != select.unique_id()) {
    equal = false;
  }
  return equal;
}

void Select::remove_last_site() {
  ASSERT(site_indices_.size() == 1, "assumes 1 particle");
  site_indices_[0].pop_back();
}

}  // namespace feasst
