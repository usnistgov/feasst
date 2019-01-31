
#include <sstream>
#include "core/include/select_position.h"

namespace feasst {

void SelectPosition::set_site_position(const int particle_index,
                                       const int site_index,
                                       const Position& position) {
  ASSERT(particle_index < static_cast<int>(particle_positions_.size()),
    "size error");
  ASSERT(site_index < static_cast<int>(site_positions_[particle_index].size()),
    "size error");
  site_positions_[particle_index][site_index] = position;
}

void SelectPosition::set_site_position(const int particle_index,
                                       const int site_index,
                                       const std::vector<double> coord) {
  Position pos;
  pos.set_vector(coord);
  set_site_position(particle_index, site_index, pos);
}

void SelectPosition::set_particle_position(const int particle_index,
                                              const Position& position) {
  ASSERT(particle_index < static_cast<int>(particle_positions_.size()),
    "size error");
  particle_positions_[particle_index] = position;
}

void SelectPosition::load_positions(const Particles& particles) {
  int pindex = 0;
  for (int particle_index : particle_indices()) {
    int sindex = 0;
    const Particle& part = particles.particle(particle_index);
    set_particle_position(pindex, part.position());
    for (int site_index : site_indices(pindex)) {
      set_site_position(pindex, sindex, part.site(site_index).position());
      site_properties_[pindex][sindex] = part.site(site_index).properties();
      ++sindex;
    }
    ++pindex;
  }
}

void SelectPosition::resize() {
  particle_positions_.resize(num_particles());
  site_positions_.resize(num_particles());
  site_properties_.resize(num_particles());
  for (int index = 0; index < num_particles(); ++index) {
    site_positions_[index].resize(site_indices(index).size());
    site_properties_[index].resize(site_indices(index).size());
  }
}

}  // namespace feasst
