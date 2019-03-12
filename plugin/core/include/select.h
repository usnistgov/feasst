
#ifndef FEASST_CORE_SELECT_H_
#define FEASST_CORE_SELECT_H_

#include <string>
#include <vector>
#include "core/include/particle_factory.h"
#include "core/include/random.h"
#include "core/include/group.h"

namespace feasst {

/**
  Select a subset of particles and sites by storing their respective indices.
  HWH Note: all input lists of indices are assumed to be sorted!
 */
class Select {
 public:
  /// Return true if nothing is selected.
  bool is_empty() const;

  /// Clear the selection.
  virtual void clear() { particle_indices_.clear(); site_indices_.clear(); }

  /// Add input selection to current.
  void add(const Select& select);

  /// Remove input selection from current.
  void remove(const Select& select);

  /// Add site by index.
  void add_site(const int particle_index, const int site_index);

  /// Add sites by index.
  void add_sites(const int particle_index,
                 const std::vector<int> site_indices);

  /// Add sites by configuration-based site index.
  void remove_sites(const int particle_index,
                    const std::vector<int> site_indices);

  /// Add particle index and site indices of that particle.
  void add_particle(const int particle_index, std::vector<int> site_indices);

//  /// Add particle by index.
//  void add_particle(const Particles& particles, const int particle_index);

  /// Add particle by index.
  void add_particle(const Particle& particle, const int particle_index);

//  /// Add last particle.
//  void add_last_particle(const Particles& particles);

  /// Remove the last particle.
  void remove_last_particle();

  /// Remove the last site.
  virtual void remove_last_site();

  /// Remove the last num sites.
  void remove_last_sites(const int num) {
    for (int index = 0; index < num; ++index) {
      remove_last_site();
    }
  }

  /// Remove the first site.
  virtual void remove_first_site();

  /// Remove the first num sites.
  void remove_first_sites(const int num) {
    for (int index = 0; index < num; ++index) {
      remove_first_site();
    }
  }

  /// Remove particle by index.
  void remove_particle(const int particle_index);

//  /// Add random particle.
//  void add_random_particle(const Particles& particles) {
//    add_particle(particles, random_particle_index(particles));
//  }

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

//  /// Return a random index of particles.
//  int random_particle_index(const Particles& particles);

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
  void check() const;

  /// Return true if the selections are equivalent.
  bool is_equivalent(const Select& select);

  /// Possible states:
  /// old -> configuration unchanged from previously accepted state
  /// add -> added new particles/sites listed in selection
  /// move -> moved selected particles but total numbers unchanged
  void set_trial_state(std::string state) { trial_state_ = state; }
  std::string trial_state() const { return trial_state_; }

//  virtual void reverse() {
//    feasst_reverse(&particle_indices_);
//    feasst_reverse(&site_indices_);
//  }

  void exclude(const Select& select) {
    if (select.num_sites() > 0) {
      excluded_ = std::make_shared<Select>();
      excluded_->add(select);
    }
  }

  const std::shared_ptr<Select> excluded() const { return excluded_; }

  virtual ~Select() {}

 protected:
  Random * random() { return &random_; }

 private:
  std::vector<int> particle_indices_;
  std::vector<std::vector<int> > site_indices_;
  Random random_;
  std::string unique_id_;
  std::string trial_state_;
  std::shared_ptr<Select> excluded_;

  // remove particle by selection index.
  void remove_particle_(const int select_index) {
    particle_indices_.erase(particle_indices_.begin() + select_index);
    site_indices_.erase(site_indices_.begin() + select_index);
  }
};

/**
  A selection based on a group.
 */
class SelectGroup : public Select {
 public:
  Group group() const { return group_; }
  void set_group(const Group group) { group_ = group; }

  virtual ~SelectGroup() {}

 private:
  Group group_;
};

}  // namespace feasst

#endif  // FEASST_CORE_SELECT_H_
