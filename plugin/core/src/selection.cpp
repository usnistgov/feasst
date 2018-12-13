
#include <sstream>
#include "core/include/selection.h"
#include "core/include/utils.h"

namespace feasst {

void Selection::add_particle(const Particles& particles,
                             const int particle_index) {
  const Particle& part = particles.particle(particle_index);
  add_particle(part, particle_index);
}

void Selection::add_particle(const Particle& particle,
                             const int particle_index) {
  std::vector<int> sites;
  for (int index = 0; index < particle.num_sites(); ++index) {
    sites.push_back(index);
  }
  particle_indices_.push_back(particle_index);
  site_indices_.push_back(sites);
}

std::string Selection::str() const {
  std::stringstream ss;
  for (int index = 0; index < num_particles(); ++index) {
    ss << "p: " << particle_indices_[index] << " s: "
       << feasst::str(site_indices_[index]) << " ";
  }
  return ss.str();
}

bool Selection::is_empty() const {
  if (num_particles() == 0) {
    return true;
  }
  return false;
}

int Selection::random_particle_index(const Particles& particles) {
  ASSERT(particles.num() > 0, "size error");
  return static_cast<int>(particles.num() * random_.uniform());
}

int Selection::num_sites() const {
  return num_elements(site_indices_);
}

Selection Selection::random_particle() {
  Selection select;
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
void Selection::reduce_to_random_particle() {
  if (num_particles() <=0) {
    return;
  }
  const int particle_index = random_.element(particle_indices_);
  const std::vector<int> sites = site_indices_[particle_index];
  clear();
  particle_indices_.push_back(particle_index);
  site_indices_.push_back(sites);
}

void Selection::add_site(const int particle_index, const int site_index) {
  int index;
  if (find_in_list(particle_index, particle_indices(), &index)) {
    site_indices_[index].push_back(site_index);
  } else {
    particle_indices_.push_back(particle_index);
    site_indices_.push_back({site_index});
  }
}

void Selection::add_sites(const int particle_index,
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

void Selection::add_last_particle(const Particles& particles) {
  add_particle(particles, particles.num() - 1);
}

void Selection::remove_last_particle() {
  particle_indices_.pop_back();
  site_indices_.pop_back();
}

void Selection::remove_particle(const int particle_index) {
  int index;
  if (find_in_list(particle_index, particle_indices_, &index)) {
    particle_indices_.erase(particle_indices_.begin() + index);
    site_indices_.erase(site_indices_.begin() + index);
  }
}

void Selection::check_size() const {
  ASSERT(particle_indices_.size() == site_indices_.size(), "size error");
}

void Selection::add_particle(const int particle_index, std::vector<int> site_indices) {
  if (site_indices.size() > 0) {
    particle_indices_.push_back(particle_index);
    site_indices_.push_back(site_indices);
  }
}

void PositionSelection::set_site_position(const int particle_index,
                                          const int site_index,
                                          const Position& position) {
  ASSERT(particle_index < static_cast<int>(particle_positions_.size()),
    "size error");
  ASSERT(site_index < static_cast<int>(site_positions_[particle_index].size()),
    "size error");
  site_positions_[particle_index][site_index] = position;
}

void PositionSelection::set_particle_position(const int particle_index,
                                              const Position& position) {
  ASSERT(particle_index < static_cast<int>(particle_positions_.size()),
    "size error");
  particle_positions_[particle_index] = position;
}

void PositionSelection::load_positions(const Particles& particles) {
  int pindex = 0;
  for (int particle_index : particle_indices()) {
    int sindex = 0;
    const Particle& part = particles.particle(particle_index);
    set_particle_position(pindex, part.position());
    for (int site_index : site_indices(pindex)) {
      set_site_position(pindex, sindex, part.site(site_index).position());
      ++sindex;
    }
    ++pindex;
  }
}

}  // namespace feasst
