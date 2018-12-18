
#ifndef FEASST_CORE_SELECT_H_
#define FEASST_CORE_SELECT_H_

#include <string>
#include <vector>
#include "core/include/particles.h"
#include "core/include/random.h"
#include "core/include/group.h"

namespace feasst {

/**
  Select a subset of particles and sites by storing their respective indices.
 */
class Select {
 public:
  /// Return true if nothing is selected.
  bool is_empty() const;

  /// Clear the selection.
  virtual void clear() { particle_indices_.clear(); site_indices_.clear(); }

  /// Add site by index.
  void add_site(const int particle_index, const int site_index);

  /// Add sites by index.
  void add_sites(const int particle_index,
                 const std::vector<int> site_indices);

  /// Add particle index and site indices of that particle.
  void add_particle(const int particle_index, std::vector<int> site_indices);

  /// Add particle by index.
  void add_particle(const Particles& particles, const int particle_index);

  /// Add particle by index.
  void add_particle(const Particle& particle, const int particle_index);

  /// Add last particle.
  void add_last_particle(const Particles& particles);

  /// Remove the last particle.
  void remove_last_particle();

  /// Remove particle by index.
  void remove_particle(const int particle_index);

  /// Add random particle.
  void add_random_particle(const Particles& particles) {
    add_particle(particles, random_particle_index(particles));
  }

  /// Return the selection of a particle chosen randomly from current
  /// selection.
  Select random_particle();

  /// Reduce existing selection to one random particle.
  void reduce_to_random_particle();

  /// Return number of selected particles.
  int num_particles() const {
    return static_cast<int>(particle_indices_.size());
  }

  /// Return number of selected sites.
  int num_sites() const;

  /// Print the selection.
  std::string str() const;

  /// Return a random index of particles.
  int random_particle_index(const Particles& particles);

  /// Return the unique identifier.
  std::string unique_id() const { return unique_id_; }

  /// Set the unique identifier.
  void set_unique_id(const std::string id) { unique_id_ = id; }

  /// Return the particle indices.
  const std::vector<int>& particle_indices() const { return particle_indices_; }

  /// Return the site indices.
  const std::vector<std::vector<int> >& site_indices() const { return site_indices_; }

  /// Return the site indices of particle
  const std::vector<int>& site_indices(const int particle_index) const {
    return site_indices_[particle_index];
  }

  int site_index(const int particle_index, const int site_index) const {
    return site_indices_[particle_index][site_index];
  }

  // Return the particle index given selection index.
  int particle_index(const int index) const { return particle_indices_[index]; }

  /// Check the size of member vectors
  void check_size() const;

 protected:
  Random * random() { return &random_; }

 private:
  std::vector<int> particle_indices_;
  std::vector<std::vector<int> > site_indices_;
  Random random_;
  std::string unique_id_;
};

/**
  A selection based on a group.
 */
class SelectGroup : public Select {
 public:
  Group group() const { return group_; }
  void set_group(const Group group) { group_ = group; }

 private:
  Group group_;
};

}  // namespace feasst

#endif  // FEASST_CORE_SELECT_H_
