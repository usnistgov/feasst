
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
  selection_.push_back(std::make_pair(particle_index, sites));
}

std::string Selection::str() const {
  std::stringstream ss;
  for (auto pair : selection_) {
    ss << "p: " << pair.first << " s: " << feasst::str(pair.second) << " ";
  }
  return ss.str();
}

bool Selection::empty() const {
  if (selection_.size() == 0) {
    return true;
  }
  return false;
}

std::vector<int> Selection::particle_indices() const {
  std::vector<int> indices;
  for (const std::pair<int, std::vector<int> >& pair : selection_) {
    indices.push_back(pair.first);
  }
  return indices;
}

int Selection::random_particle_index(const Particles& particles) {
  ASSERT(particles.num() > 0, "size error");
  return static_cast<int>(particles.num() * random_.uniform());
}

int Selection::num_sites() const {
  int num = 0;
  for (const std::pair<int, std::vector<int> >& pair : selection_) {
    num += static_cast<int>(pair.second.size());
  }
  return num;
}

void Selection::reduce_to_random_particle() {
  if (num_particles() <=0) {
    return;
  }
  selection_[0] = random_.element(selection_);
}

void Selection::add_site(const int particle_index, const int site_index) {
  int index;
  if (find_in_list(particle_index, particle_indices(), &index)) {
    ERROR("HWH implement this");
  } else {
    std::vector<int> sites = {site_index};
    selection_.push_back(std::make_pair(particle_index, sites));
  }
}

}  // namespace feasst
