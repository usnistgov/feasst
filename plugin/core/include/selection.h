
#ifndef FEASST_CORE_SELECTION_H_
#define FEASST_CORE_SELECTION_H_

#include <string>
#include <vector>
#include <utility>  // pair
#include "core/include/particles.h"
#include "core/include/random.h"

namespace feasst {

/**
  Select a subset of particles and sites by storing their respective indices.
 */
class Selection {
 public:
  /// Return true if nothing is selected.
  bool empty() const;

  /// Clear the selection.
  void clear() { selection_.clear(); }

  /// Add site by index.
  void add_site(const int particle_index, const int site_index);

  /// Add particle by index.
  void add_particle(const Particles& particles, const int particle_index);

  /// Add particle by index.
  void add_particle(const Particle& particle, const int particle_index);

  /// Add random particle.
  void add_random_particle(const Particles& particles) {
    add_particle(particles, random_particle_index(particles));
  }

  /// Reduce existing selection to one random particle.
  void reduce_to_random_particle();

  /// Return number of selected particles.
  int num_particles() const { return static_cast<int>(selection_.size()); }

  /// Return number of selected sites.
  int num_sites() const;

  /// Print the selection.
  std::string str() const;

  /// Return the particle indices.
  std::vector<int> particle_indices() const;

  /// Return a random index of particles.
  int random_particle_index(const Particles& particles);

  /// Return the unique identifier.
  std::string unique_id() const { return unique_id_; }

  /// Set the unique identifier.
  void set_unique_id(const std::string id) { unique_id_ = id; }

  /*
    Return the selection as a list of the following pair of quantities:

      1. particle index.
      2. list of site indices of that particle.

    Note HWH: discourage this usage? how to best use selection?
   */
  const std::vector<std::pair<int, std::vector<int> > >& selection() const {
    return selection_; }

  // Return the particle index given selection index.
  // Note HWH: discourage this usage?
  int particle_index(const int index) const { return selection_[index].first; }

 private:
  std::vector<std::pair<int, std::vector<int> > > selection_;
  Random random_;
  std::string unique_id_;
};

}  // namespace feasst

#endif  // FEASST_CORE_SELECTION_H_
