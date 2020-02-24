
#include <sstream>
#include "configuration/include/select_position.h"

namespace feasst {

SelectPosition::SelectPosition(const Select& select,
  const ParticleFactory& particles)
  : Select(select) {
  resize();
  load_positions(particles);
}

void SelectPosition::set_site_position(const int particle_index,
                                       const int site_index,
                                       const Position& position) {
  ASSERT(particle_index < static_cast<int>(particle_positions_.size()),
    "particle_index(" << particle_index << ") is larger than position size: "
    << particle_positions_.size());
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

void SelectPosition::add_to_site_position(const int particle_index,
                                          const int site_index,
                                          const Position& position) {
  ASSERT(particle_index < static_cast<int>(particle_positions_.size()),
    "particle_index(" << particle_index << ") is larger than position size: "
    << particle_positions_.size());
  ASSERT(site_index < static_cast<int>(site_positions_[particle_index].size()),
    "size error");
  site_positions_[particle_index][site_index].add(position);
}

void SelectPosition::add_to_particle_position(const int particle_index,
                                              const Position& position) {
  ASSERT(particle_index < static_cast<int>(particle_positions_.size()),
    "particle_index(" << particle_index << ") is larger than position size: "
    << particle_positions_.size());
  particle_positions_[particle_index].add(position);
}

void SelectPosition::set_site_properties(
    const int particle_index,
    const int site_index,
    const Properties& properties) {
  site_properties_[particle_index][site_index] = properties;
}

void SelectPosition::set_particle_position(const int particle_index,
                                              const Position& position) {
  ASSERT(particle_index < static_cast<int>(particle_positions_.size()),
    "size error");
  particle_positions_[particle_index] = position;
}

void SelectPosition::load_position(const int pindex,
    const Particle& particle) {
  set_particle_position(pindex, particle.position());
  int sindex = 0;
  for (int site_index : site_indices(pindex)) {
    set_site_position(pindex, sindex, particle.site(site_index).position());
    set_site_properties(pindex, sindex, particle.site(site_index).properties());
    ++sindex;
  }
}

void SelectPosition::load_positions(const ParticleFactory& particles) {
  int pindex = 0;
  for (int particle_index : particle_indices()) {
    load_position(pindex, particles.particle(particle_index));
    ++pindex;
  }
}

void SelectPosition::load_positions_of_last(const Particle& particle,
    const Position& frame_of_reference) {
  particle_positions_.push_back(particle.position());
  particle_positions_.back().add(frame_of_reference);
  const int pindex = num_particles() - 1;
  std::vector<Position> site_posit;
  std::vector<Properties> site_props;
  for (int sindex : site_indices(pindex)) {
    site_posit.push_back(particle.site(sindex).position());
    site_posit.back().add(frame_of_reference);
    site_props.push_back(particle.site(sindex).properties());
  }
  site_positions_.push_back(site_posit);
  site_properties_.push_back(site_props);
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

void SelectPosition::remove_last_site() {
  Select::remove_last_site();
  site_positions_[0].pop_back();
  site_properties_[0].pop_back();
}

void SelectPosition::remove_first_site() {
  Select::remove_first_site();
  site_positions_[0].erase(site_positions_[0].begin());
  site_properties_[0].erase(site_properties_[0].begin());
}

void SelectPosition::clear() {
  Select::clear();
  clear_();
}

void SelectPosition::add_site(const int particle_index, const int site_index) {
  Select::add_site(particle_index, site_index);
  resize();
}

void SelectPosition::clear_() {
  particle_positions_.clear();
  site_positions_.clear();
  site_properties_.clear();
}

Position SelectPosition::geometric_center(const int particle_index) const {
  Position center(particle_positions()[0].dimension());
  // consider all particles if particle index is not provided
  if (particle_index == -1) {
    for (int sp = 0; sp < num_particles(); ++sp) {
      for (int ss = 0; ss < num_sites(sp); ++ss) {
        center.add(site_positions()[sp][ss]);
      }
    }
    center.divide(num_sites());
  } else {
    const int sp = particle_index;
    for (int ss = 0; ss < num_sites(sp); ++ss) {
      center.add(site_positions()[sp][ss]);
    }
    center.divide(num_sites(particle_index));
  }
  return center;
}

SelectPosition::SelectPosition(const int particle_index,
    const Particle& particle) {
  add_particle(particle, particle_index);
  resize();
  load_position(0, particle);
}

void SelectPosition::serialize(std::ostream& sstr) const {
  feasst_serialize_version(800, sstr);
  feasst_serialize_fstobj(particle_positions_, sstr);
  feasst_serialize_fstobj(site_positions_, sstr);
  feasst_serialize_fstobj(site_properties_, sstr);
}

SelectPosition::SelectPosition(std::istream& sstr) {
  const int version = feasst_deserialize_version(sstr);
  ASSERT(version == 800, "version");
  feasst_deserialize_fstobj(&particle_positions_, sstr);
  feasst_deserialize_fstobj(&site_positions_, sstr);
  feasst_deserialize_fstobj(&site_properties_, sstr);
}
}  // namespace feasst
