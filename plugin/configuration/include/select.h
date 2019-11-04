
#ifndef FEASST_CONFIGURATION_SELECT_H_
#define FEASST_CONFIGURATION_SELECT_H_

#include <string>
#include <vector>
#include "configuration/include/particle_factory.h"
#include "math/include/random.h"
#include "configuration/include/group.h"

namespace feasst {

/**
  Select a subset of particles and sites by storing their respective indices.
  For optimization, all input lists of indices are assumed to be sorted.
 */
class Select {
 public:
  Select() {
    set_trial_state();
  }

  /// Return true if nothing is selected.
  bool is_empty() const;

  /// Clear the selection.
  virtual void clear() { particle_indices_.clear(); site_indices_.clear(); }

  /// Add input selection to current.
  void add(const Select& select);

  /// Remove input selection from current.
  void remove(const Select& select);

  /// Add site by index.
  virtual void add_site(const int particle_index, const int site_index);

  /// Set site by index
  void set_site(const int particle_index,
      const int site_index,
      const int index) {
    site_indices_[particle_index][site_index] = index; }

  /// Set particle by index
  void set_particle(const int particle_index, const int index) {
    particle_indices_[particle_index] = index; }

  /// Add sites by index.
  void add_sites(const int particle_index,
                 const std::vector<int> site_indices);

  /// Remove sites by configuration-based site index.
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
  void remove_last_sites(const int num);

  /// Remove the first site.
  virtual void remove_first_site();

  /// Remove the first num sites.
  void remove_first_sites(const int num);

  /// Remove particle by index.
  void remove_particle(const int particle_index);

  /// Return the selection of a particle chosen randomly from current
  /// selection.
  Select random_particle(Random * random);

  /// Return number of selected particles.
  int num_particles() const {
    return static_cast<int>(particle_indices_.size());
  }

  /// Return number of selected sites.
  int num_sites(
    /// Return numbers of sites in particle. If -1 (default), in all particles.
    const int particle_index = -1) const;

  /// Print the selection.
  std::string str() const;

  /// Return the particle indices.
  const std::vector<int>& particle_indices() const { return particle_indices_; }

  /// Return the site indices.
  const std::vector<std::vector<int> >& site_indices() const {
    return site_indices_; }

  /// Return the site indices of particle
  const std::vector<int>& site_indices(const int particle_index) const {
    return site_indices_[particle_index];
  }

  int site_index(const int particle_index, const int site_index) const {
    return site_indices_[particle_index][site_index]; }

  // Return the particle index given selection index.
  int particle_index(const int index) const { return particle_indices_[index]; }

  /// Check the size of member vectors
  void check() const;

  /// Return true if the selections are equivalent.
  bool is_equal(const Select& select) const;

  /// Possible states:
  /// 0 -> old -> configuration unchanged from previously accepted state
  /// 1 -> move -> moved selected particles but total numbers unchanged
  /// 2 -> remove -> remove existing particles/sites listed in selection
  /// 3 -> add -> added new particles/sites listed in selection
  int trial_state() const { return trial_state_; }

  /// Set the trial state.
  void set_trial_state(const int state = -1) { trial_state_ = state; }

//  virtual void reverse() {
//    feasst_reverse(&particle_indices_);
//    feasst_reverse(&site_indices_);
//  }

  /// Excluded site list.
  void exclude(const Select& select);

  /// Return excluded list.
  const std::shared_ptr<Select> excluded() const { return excluded_; }

  /// Sites to become bonded.
  void set_new_bond(const Select& select);

  /// Return sites to become bonded.
  const std::shared_ptr<Select> new_bond() const { return new_bond_; }

  /// Sites which are bonded.
  void set_old_bond(const Select& select);

  /// Return sites which are bonded.
  const std::shared_ptr<Select> old_bond() const { return old_bond_; }

  /// Reset excluded and bonded sites.
  void reset_excluded_and_bond();

  /// Replace current indices with those given. Return true if replace is done
  /// quickly due to match in existing size.
  bool replace_indices(const int particle_index,
    const std::vector<int>& site_indices);

  /// Return true if the selection contains a particle in self.
  bool is_overlap(const Select& select) const;

  virtual void serialize(std::ostream& ostr) const;
  Select(std::istream& istr);
  virtual ~Select() {}

 private:
  std::vector<int> particle_indices_;
  std::vector<std::vector<int> > site_indices_;
  int trial_state_;
  std::shared_ptr<Select> excluded_;
  std::shared_ptr<Select> new_bond_;
  std::shared_ptr<Select> old_bond_;

  // remove particle by selection index.
  void remove_particle_(const int select_index);
};

/**
  A selection based on a group.
 */
class SelectGroup : public Select {
 public:
  SelectGroup() {}
  const Group& group() const { return group_; }
  void set_group(const Group group) { group_ = group; }

  void serialize(std::ostream& ostr) const override;
  SelectGroup(std::istream& istr);
  virtual ~SelectGroup() {}

 private:
  Group group_;
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_SELECT_H_
