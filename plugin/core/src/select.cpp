
#include <sstream>
#include "core/include/select.h"
#include "core/include/utils.h"

namespace feasst {

void Select::add_particle(const Particles& particles,
                             const int particle_index) {
  const Particle& part = particles.particle(particle_index);
  add_particle(part, particle_index);
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
    ss << "p: " << particle_indices_[index] << " s: "
       << feasst::str(site_indices_[index]) << " ";
  }
  return ss.str();
}

bool Select::is_empty() const {
  if (num_particles() == 0) {
    return true;
  }
  return false;
}

int Select::random_particle_index(const Particles& particles) {
  ASSERT(particles.num() > 0, "size error");
  return static_cast<int>(particles.num() * random_.uniform());
}

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
    std::vector<int> * sind = &site_indices_[index];
    sind->insert(sind->end(), site_indices.begin(), site_indices.end());
  } else {
    particle_indices_.push_back(particle_index);
    site_indices_.push_back(site_indices);
  }
}

void Select::add_last_particle(const Particles& particles) {
  add_particle(particles, particles.num() - 1);
}

void Select::remove_last_particle() {
  particle_indices_.pop_back();
  site_indices_.pop_back();
}

void Select::remove_particle(const int particle_index) {
  int index;
  if (find_in_list(particle_index, particle_indices_, &index)) {
    particle_indices_.erase(particle_indices_.begin() + index);
    site_indices_.erase(site_indices_.begin() + index);
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

}  // namespace feasst
