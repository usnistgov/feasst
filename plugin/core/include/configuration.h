
#ifndef FEASST_CORE_CONFIGURATION_H_
#define FEASST_CORE_CONFIGURATION_H_

#include <string>
#include <vector>
#include "core/include/domain.h"
#include "core/include/particles.h"
#include "core/include/select.h"
#include "core/include/select_position.h"

namespace feasst {

/**
  A configuration contains both the particles and the spatial domain/boundaries.

  For the particles, this includes both the pool of particles which may
  exist (referred to as particle types) in addition to the particles which
  physically exist (referred to as particles).

  Before particles are physically added to the system, all particle types which
  may exist must be defined.
  Once particles are added to the system, new particle types cannot be added.
  Similarly, site types are also defined at the same time as particle types.

  Groups of different particle/site types and other metrics may be defined.
  These groups then define a selection which can be used to distinguish subsets
  of the configuration (e.g., types of particles).
  This selection may be further reduced to single particles.
  These selections are then used to modify a subset of the configuration (e.g.,
  removal and displacement) of a selection of particles/sites.

  Finally, the spatial domain contains periodic boundaries and cells.
 */
class Configuration {
 public:
  Configuration();

  /** @name Typing
    Types of sites and particles. */
  //@{

  /// Add a particle type that may exist by LMP file (see FileLMP).
  void add_particle_type(const char* file_name);

  /// Return the number of particle types.
  int num_particle_types() const { return particle_types_.num(); }

  /// Return the number of site types.
  int num_site_types() const { return unique_types_.num_sites(); }

  /// Return the particle types.
  const Particles& particle_types() const { return particle_types_; }

  /// Modify model parameter of a given site type and name to value.
  void set_model_param(const char* name,
                       const int site_type,
                       const double value) {
    unique_types_.set_model_param(name, site_type, value);
  }

  /// Add model parameter of a given name to value.
  void add_model_param(const std::string name,
                       const double value) {
    unique_types_.add_model_param(name, value);
  }

  /// Return the particle associated with the type.
  const Particle& particle_type(const int type) const {
    return particle_types_.particle(type);
  }

  /// Return the file name used to initialize the particle types.
  std::string type_to_file_name(const int type) const {
    return type_to_file_[type];
  }

  /// Return the unique types. Only unique sites and bonds are included.
  /// Thus, the site index is the same as the numeric value for the site type.
  /// And the same for bonds.
  /// This serves as a container for properties based on site or bond type.
  const Particles& unique_types() const { return unique_types_; }

  /// Return the unique type by individual particle.
  const Particle& unique_type(const int type) const {
    return unique_types_.particle(type);
  }

  //@}
  /** @name Groups
    Groups of sites and particles
   */
  //@{

  /// Add a group (after types are defined but before particles are added).
  void add(Group group,
    /// Optionally provide a name. If no name is provided, the name is assigned
    /// to be the numerical indices of the order of groups added.
    std::string name = "");

  /// Return the number of group selections.
  int num_groups() const { return static_cast<int>(group_selects_.size()); }

  /// Return the index of the group based on particle types.
  int particle_type_to_group(const int particle_type);

  //@}
  /** @name Particles
    Physically existing sites and particles
   */
  //@{

  /// Add particles of a given type.
  void add_particle(const int type);

  /// Return particle by index. Note this index is contiguous from values
  /// 0 to num_particles -1, unlike the selection indices (due to ghost
  /// particles).
  /// Note that this method can be quite slow because it doesn't require
  /// knowledge of ghost particles.
  const Particle particle(const int index,
    /// Provide a group index to consider only a subset of the configuration.
    /// By default, a value of zero is for the entire configuration.
    const int group = 0) const;

  /// Return the number of particles.
  int num_particles(
    /// Provide a group index to consider only a subset of the configuration.
    /// By default, a value of zero is for the entire configuration.
    const int group = 0) const;

  /// Return the number of sites.
  int num_sites() const { return particles_.num_sites(); }

  //@}
  /** @name Modifications
    Modifications to a configuration (e.g., moving, adding or deleting
    particles/sites.
    A subset of the configuration is defined by a Select.
    Avoid using a Select or SelectPosition object directly.
    Instead, use a SelectList, which can be input into the following.
   */
  //@{

  /// Load coordinates by per-site vector containing per-dimension vector.
  /// Requires coordinates for all sites and dimensions.
  void update_positions(const std::vector<std::vector<double> > coords);

  /// Update the positions from a selection.
  void update_positions(const SelectPosition& select);

  /// Displace selected particle(s).
  void displace_particles(const Select& selection,
                          const Position &displacement);

  /// Same as above except for only one particle that is selected.
  void displace_particle(const Select& selection, const Position &displacement);

  /// Displace the selection. No periodic boundary conditions applied.
  void displace(const Select& selection, const Position &displacement);

  /// Replace positions of particle by selection.
  void replace_position(const Select& select, const Particle& replacement);

  /// Remove particle(s) in selection.
  void remove_particles(const Select& selection);

  /// Same as above except for only one particle that is selected.
  void remove_particle(const Select& selection);

  /// Revive the selection which was deleted.
  void revive(const SelectPosition& selection);

  //@{
  /** @name Domain
    A configuration's domain includes periodic boundaries and cells.
   */
  //${

  /// Return the domain of the configuration.
  const Domain& domain() const { return domain_; }

  /// Set the domain side lengths.
  // HWH consider scaling particles as well
  void set_side_length(const Position& sides) {
//    ASSERT(num_particles() == 0, "domain should only be scaled with particles");
    domain_.set_side_length(sides);
  }

  /// Set the domain.
  // HWH depreciate
  void set_domain(const Domain domain) { domain_ = domain; }

  /// Return the dimensionality of space.
  int dimension() const { return domain().dimension(); }

  /// Initialize the cells according to the minimum side length.
  void init_cells(const double min_length,
    /// By default, cells are applied to all particles and sites.
    /// Set the group index to consider only a subset.
    const int group_index = 0);

  //@}
  /*
    Functions which require knowledge of ghost particles and thus not for
    typical users
   */

  // Return selection of all particles and sites in the configuration.
  const Select& selection_of_all() const { return group_selects_[0]; }

  // Return the particles.
  // Warning: typically not for users because it may include ghost particles.
  const Particles& particles() const { return particles_; }

  // Return particle by index provided in selection.
  // Warning: typically not for users because it may include ghost particles.
  const Particle& select_particle(const int index) const {
    return particles_.particle(index); }

  // Return the selection-based index (includes ghosts) of the last particle added.
  int newest_particle_index() const { return newest_particle_index_; }
  const Particle& newest_particle() const {
    return select_particle(newest_particle_index_); }

  // Return the group-based selections.
  const std::vector<SelectGroup>& group_selects() const {
    return group_selects_; }

  // Return the group-based selections by index.
  const SelectGroup& group_select(const int index) const {
    return group_selects_[index]; }

  // Add the property to a site in a particle.
  void add_site_property(const std::string name,
      const double value,
      const int particle_index,
      const int site_index) {
    particles_.add_site_property(name, value, particle_index, site_index);
  }

  // Set the property of a site in a particle.
  void set_site_property(const std::string name,
      const double value,
      const int particle_index,
      const int site_index) {
    particles_.set_site_property(name, value, particle_index, site_index);
  }
  void set_site_property(const int index,
      const double value,
      const int particle_index,
      const int site_index) {
    particles_.set_site_property(index, value, particle_index, site_index);
  }

  /* Checks and hacky additions */

  /// Check consistency of dimensions and lists.
  void check_size() const;

 private:
  Particles particle_types_;
  Particles unique_types_;
  Particles particles_;
  Domain domain_;
  int newest_particle_index_;

  /// Selects based on groups that are continuously updated.
  // HWH currently only updated when adding and removing particles
  // HWH but at some point it should check for positional changes
  // HWH if groups are defined based on positions.
  std::vector<SelectGroup> group_selects_;

  /// Unique identifier for the collection of particle indices.
  std::string unique_indices_;
  void reset_unique_indices_();
  Random random_;

  /// Add particle.
  void add_(const Particle particle);

  /// Remove particle by index.
  void remove_particle_(const int index);

  /// Displace particle by index.
  void displace_particle_(const int index,
                          const Position &displacement);

  /// Displace site by index.
  void displace_site_(const int particle_index,
                      const int site_index,
                      const Position &displacement);

  /// Replace positions of particle by index.
  void replace_position_(const int index,
                         const Particle& replacement);

  /// Replace position of site in particle.
  void replace_position_(const int particle_index,
                         const int site_index,
                         const Position& replacement);

  /// Replace position of particle but not site.
  void replace_position_(const int particle_index,
                         const Position& replacement);

  /// Replace properties of site in particle.
  void replace_properties_(const int particle_index,
                           const int site_index,
                           const Properties& properties) {
    particles_.replace_properties(particle_index, site_index, properties); }

  /// Update position trackers of a particle (e.g., cell, neighbor, etc).
  void position_tracker_(const int particle_index, const int site_index);

  /// Update position trackers of all sites in a particle.
  void position_tracker_(const int particle_index);

  /// Update position trackers of all particles.
  void position_tracker_();

  /// Add particle to selection.
  void add_to_selection_(const int particle_index,
                         SelectGroup * select) const;

  /// Initialize selection based on groups
  void init_selection_(SelectGroup * group_select) const;

  /// Remember groups based on types.
  std::vector<int> group_store_particle_type_,
                   group_store_group_index_;

  /// HWH depreciate one of these.
  void check_id_(const Select& select) const;
  void check_id_(const Select* select) const;
  void check_id_(const std::string id) const;

  // Ghost particles allow quick addition and deletion of particles for use in
  // the grand canonical ensemble.
  // ghost particles types are related to real particle types as follows
  // type_ghost = -type - 1
  // ghosts are removed from selections and can be brought back by adding.
  std::vector<SelectGroup> ghosts_;

  /// Return the number of ghost particles.
  int num_ghosts_() const;

  const Particle& particle_(const int index) {
    return particles_.particle(index);
  }

  /// Store the files used to initialize particle types.
  std::vector<std::string> type_to_file_;
};

}  // namespace feasst

#endif  // FEASST_CORE_CONFIGURATION_H_
