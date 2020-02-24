
#include <sstream>
#include <algorithm>
#include "configuration/include/select.h"
#include "utils/include/utils.h"
#include "math/include/utils_math.h"
#include "utils/include/utils_io.h"
#include "utils/include/debug.h"

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
  TRACE("removing ");
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
  add_particle(particle_index, sites);
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

int Select::num_sites(const int particle_index) const {
  if (particle_index == -1) {
    return num_elements(site_indices_);
  }
  ASSERT(particle_index < static_cast<int>(site_indices_.size()),
    "particle_index:" << particle_index << " is too large for selection "
    << " with " << site_indices_.size() << " particles");
  return static_cast<int>(site_indices_[particle_index].size());
}

Select Select::random_particle(Random * random) {
  Select select;
  if (num_particles() > 0) {
    int select_index;
    const int particle_index = random->const_element(particle_indices_,
                                                       &select_index);
    DEBUG("num " << num_particles());
    DEBUG("index " << particle_index);
    const std::vector<int> sites = site_indices_[select_index];
    select.add_particle(particle_index, sites);
  }
  return select;
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
    site_indices_[index] = feasst_union(site_indices_[index], site_indices);
  } else {
    particle_indices_.push_back(particle_index);
    site_indices_.push_back(site_indices);
  }
}

void Select::remove_sites(const int particle_index,
                          const std::vector<int> site_indices) {
  TRACE("removing particle_index " << particle_index);
  TRACE("removing site_indices " << feasst_str(site_indices));
  int index;
  if (find_in_list(particle_index, particle_indices(), &index)) {
    site_indices_[index] = feasst_difference(site_indices_[index], site_indices);
    // completely remove particle if no more sites remain.
    if (site_indices_[index].size() == 0) {
      remove_particle_(index);
    }
  }
}

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

void Select::check() const {
  ASSERT(particle_indices_.size() == site_indices_.size(), "size error");
  for (auto sites : site_indices_) {
    ASSERT(std::is_sorted(sites.begin(), sites.end()), "must be sorted");
  }
}

void Select::add_particle(const int particle_index, std::vector<int> site_indices) {
  if (site_indices.size() > 0) {
    particle_indices_.push_back(particle_index);
    site_indices_.push_back(site_indices);
  }
}

bool Select::is_equal(const Select& select) const {
  if (particle_indices_ != select.particle_indices()) {
    DEBUG("particle indices are not equivalent: " <<
      feasst_str(particle_indices_) << " vs " <<
      feasst_str(select.particle_indices()));
    return false;
  }
  if (site_indices_ != select.site_indices()) {
    DEBUG("site indices are not equivalent: " <<
      feasst_str(site_indices_) <<
      feasst_str(select.site_indices()));
    return false;
  }
  return true;
}

void Select::remove_last_site() {
  ASSERT(site_indices_.size() == 1, "assumes 1 particle");
  site_indices_[0].pop_back();
}

void Select::remove_last_sites(const int num) {
  for (int index = 0; index < num; ++index) {
    remove_last_site();
  }
}

void Select::remove_first_site() {
  ASSERT(site_indices_.size() == 1, "assumes 1 particle");
  site_indices_[0].erase(site_indices_[0].begin());
}

bool Select::is_overlap(const Select& select) const {
  for (const int part : select.particle_indices()) {
    if (find_in_list(part, particle_indices())) {
      return true;
    }
  }
  return false;
}

void Select::remove_first_sites(const int num) {
  for (int index = 0; index < num; ++index) {
    remove_first_site();
  }
}

void Select::exclude(const Select& select) {
  if (select.num_sites() > 0) {
    excluded_ = std::make_shared<Select>();
    excluded_->add(select);
  }
}

void Select::set_new_bond(const Select& select) {
  if (select.num_sites() > 0) {
    new_bond_ = std::make_shared<Select>();
    new_bond_->add(select);
  }
}

void Select::set_old_bond(const Select& select) {
  if (select.num_sites() > 0) {
    old_bond_ = std::make_shared<Select>();
    old_bond_->add(select);
  }
}

void Select::reset_excluded_and_bond() {
  excluded_.reset();
  new_bond_.reset();
  old_bond_.reset();
}

bool Select::replace_indices(const int particle_index,
    const std::vector<int>& site_indices) {
  if (static_cast<int>(particle_indices_.size()) == 1 and
      site_indices_.size() == site_indices.size()) {
    particle_indices_[0] = particle_index;
    site_indices_[0] = site_indices;
    return true;
  }
  clear();
  add_particle(particle_index, site_indices);
  return false;
}

void Select::remove_particle_(const int select_index) {
  particle_indices_.erase(particle_indices_.begin() + select_index);
  site_indices_.erase(site_indices_.begin() + select_index);
}

void Select::serialize(std::ostream& sstr) const {
  feasst_serialize_version(183, sstr);
  feasst_serialize(particle_indices_, sstr);
  feasst_serialize(site_indices_, sstr);
  feasst_serialize(trial_state_, sstr);
  feasst_serialize(excluded_, sstr);
  feasst_serialize(new_bond_, sstr);
  feasst_serialize(old_bond_, sstr);
}

Select::Select(std::istream& sstr) {
  const int version = feasst_deserialize_version(sstr);
  ASSERT(version == 183, "version mismatch: " << version);
  feasst_deserialize(&particle_indices_, sstr);
  feasst_deserialize(&site_indices_, sstr);
  feasst_deserialize(&trial_state_, sstr);
  feasst_deserialize(excluded_, sstr);
  feasst_deserialize(new_bond_, sstr);
  feasst_deserialize(old_bond_, sstr);
}

void SelectGroup::serialize(std::ostream& sstr) const {
  Select::serialize(sstr);
  feasst_serialize_version(1, sstr);
  group_.serialize(sstr);
}

SelectGroup::SelectGroup(std::istream& sstr) : Select(sstr) {
  feasst_deserialize_version(sstr);
  group_ = Group(sstr);
}

}  // namespace feasst
