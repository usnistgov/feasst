
#ifndef FEASST_CONFIGURATION_SELECT_H_
#define FEASST_CONFIGURATION_SELECT_H_

#include <string>
#include <vector>
#include <memory>
#include "configuration/include/group.h"
#include "configuration/include/particle_factory.h"

namespace feasst {

class Random;

// HWH rename, many particle_index are actually select_index
// HWH consider refactoring selection
// HWH could use std:pair to couple particle_index with site_indices
// HWH or use a std::set or std::multiset to enforce order?
/**
  Select a subset of particles and sites by storing their respective indices.
  For optimization, all input lists of indices are assumed to be sorted.
 */
class Select {
 public:
  //@{
  /** @name Construction
   */

  Select() {
    set_trial_state();
  }

  //@}
  /** @name Indices
    Add, remove or access the indices of sites and particles.
   */
  //@{

  /// Return true if nothing is selected.
  bool is_empty() const;

  /// Clear the selection.
  void clear();

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
  void add_sites(
    const int particle_index,
    const std::vector<int> site_indices);

  /// Remove sites by configuration-based particle and site index.
  void remove_sites(const int particle_index,
                    const std::vector<int> site_indices);

  /// Add particle index and site indices of that particle.
  void add_particle(const int particle_index,
    std::vector<int> site_indices,
    /// optionally prevent duplicates (slower)
    const bool prevent_duplicate = false);

  /// Same as above, but use Particle to get site_indices.
  void add_particle(const Particle& particle,
    const int particle_index,
    const bool prevent_duplicate = false);

  /// Return number of selected particles.
  int num_particles() const {
    return static_cast<int>(particle_indices_.size());
  }

  /// Return number of selected sites.
  int num_sites(
    /// Return numbers of sites in particle. If -1 (default), in all particles.
    const int particle_index = -1) const;

//  /// Add last particle.
//  void add_last_particle(const Particles& particles);

  /// Add input selection to current.
  void add(const Select& select);

  /// Remove input selection from current.
  void remove(const Select& select);

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

  /// Return the particle indices.
  const std::vector<int>& particle_indices() const { return particle_indices_; }

  /// Return the site indices.
  const std::vector<std::vector<int> >& site_indices() const {
    return site_indices_; }

  /// Return the site indices of particle
  /// Note that particle_index is the index of selection, not particle index
  ///   from configuration.
  const std::vector<int>& site_indices(const int particle_index) const {
    return site_indices_[particle_index];
  }

  int site_index(const int particle_index, const int site_index) const {
    return site_indices_[particle_index][site_index]; }

  // Return the particle index given selection index.
  int particle_index(const int index) const { return particle_indices_[index]; }

  /// Replace current indices with those given. Return true if replace is done
  /// quickly due to match in existing size.
  bool replace_indices(const int particle_index,
    const std::vector<int>& site_indices);

  /// Return true if the selection contains a particle in self.
  bool is_overlap(const Select& select) const;

  //@}
  /** @name Group
    Group information.
   */
  //@{

  /// Return true if group is defined.
  bool is_group_empty() const { if (group_) { return false; } return true; }

  /// Set the group.
  void set_group(std::shared_ptr<Group> group) { group_ = group; }

  /// Return the group.
  const Group& group() const;

  //@}
  /** @name Positions
    Positions and properties.
   */
  //@{

  /// Construct with positions.
  Select(const Select& select, const ParticleFactory& particles);

  /// Construct with positions.
  Select(const int particle_index, const Particle& particle);

  /// Return true if select has positions.
  bool has_positions() const;

  // HWH move this to private to perform upon load_positions
  void resize_positions();

  /// Return the site positions.
  const std::vector<std::vector<Position> >& site_positions() const {
    return site_positions_; }

  /// Return the site positions.
  const std::vector<std::vector<Properties> >& site_properties() const {
    return site_properties_; }

  /// Set the position of a site by particle and site index.
  /// Note that these indices are based on selection, not configuration.
  void set_site_position(const int particle_index,
                         const int site_index,
                         const Position& position);

  /// Same as above except vector position is accepted.
  void set_site_position(const int particle_index,
                         const int site_index,
                         const std::vector<double> coord);

  /// Add to the position of a site by particle and site index.
  /// Note that these indices are based on selection, not configuration.
  void add_to_site_position(const int particle_index,
                            const int site_index,
                            const Position& position);

  /// Set the property of a site by particle and site index.
  /// Note that these indices are based on selection, not configuration.
  void set_site_properties(const int particle_index,
                           const int site_index,
                           const Properties& properties);

  /// Load the positions of a particle with existing selection indices.
  void load_position(const int pindex, const Particle& particle);

  /// Load the positions from the existing selection indices.
  void load_positions(const ParticleFactory& particles);

  /// Load the positions and properties of the last particle added.
  void load_positions_of_last(const Particle& particle,
                              /// shift the positions by the frame of reference
                              const Position& frame_of_reference);

  /// Return the geometric center of the selection.
  Position geometric_center(
    /// Consider only one particle, or all particles (-1).
    const int particle_index = -1) const;

  // optimized access to site positions.
  Position * get_site_position(const int particle_index, const int site_index) {
    return &site_positions_[particle_index][site_index]; }

  //@}
  /** @name Trials
    Advanced trial information: state, exclude, include.
   */
  //@{

  /// Possible states:
  /// 0 -> old -> configuration unchanged from previously accepted state
  /// 1 -> move -> moved selected particles but total numbers unchanged
  /// 2 -> remove -> remove existing particles/sites listed in selection
  /// 3 -> add -> added new particles/sites listed in selection
  int trial_state() const { return trial_state_; }

  /// Set the trial state.
  void set_trial_state(const int state = -1);

  /// Excluded site list.
  void exclude(const Select& select);

  /// Return excluded list.
  const std::shared_ptr<Select> excluded() const { return excluded_; }

  /// Sites to become bonded.
  void set_new_bond(const Select& select);

  /// Return sites to become bonded.
  const Select * new_bond() const { return new_bond_.get(); }
  //const Select& new_bond() const { return *new_bond_; }

  /// Sites which are bonded.
  void set_old_bond(const Select& select);

  /// Return sites which are bonded.
  const Select * old_bond() const { return old_bond_.get(); }
  // const Select& old_bond() const { return *old_bond_; }

  /// Reset excluded and bonded sites.
  void reset_excluded_and_bond();

  //@}
  /** @name Checks
    Consistency checks and tests.
   */
  //@{

  /// Check the size of member vectors
  void check() const;

  /// Return true if the selections are equivalent.
  bool is_equal(const Select& select) const;

  //@}

  /// Return the selection of a particle chosen randomly from current
  /// selection.
  Select random_particle(Random * random);

  /// Print the selection.
  std::string str() const;

  virtual void serialize(std::ostream& ostr) const;
  explicit Select(std::istream& istr);
  virtual ~Select() {}

 private:
  std::vector<int> particle_indices_;
  std::vector<std::vector<int> > site_indices_;
  int trial_state_;
  std::shared_ptr<Select> excluded_;
  std::shared_ptr<Select> new_bond_;
  std::shared_ptr<Select> old_bond_;
  std::shared_ptr<Group> group_;
  std::vector<std::vector<Position> > site_positions_;
  std::vector<std::vector<Properties> > site_properties_;

  // remove particle by selection index.
  void remove_particle_(const int select_index);
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_SELECT_H_
