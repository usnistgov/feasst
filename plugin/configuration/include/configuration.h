
#ifndef FEASST_CONFIGURATION_CONFIGURATION_H_
#define FEASST_CONFIGURATION_CONFIGURATION_H_

#include <memory>
#include <string>
#include <vector>
#include "utils/include/arguments.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/select.h"

namespace feasst {

class Domain;

/**
  A configuration contains both the particles and the spatial domain/boundaries.

  For the particles, this includes both the pool of particles which may
  exist (referred to as particle types) in addition to the particles which
  physically exist (referred to as particles).

  Before particles are physically added to the system, all particle types which
  may exist must be defined.
  Once particles are added to the system, new particle types cannot be added.
  Similarly, site and bond types are also defined at the same time as particle
  types and are referred to as unique types.
  Unique types are stored in essentially the same fashion as particle types,
  except they have been stripped of their non-unique sites and bonds.
  Two different particles cannot share a site type.

  Groups of different particle/site types and other metrics may be defined.
  These groups then define a selection which can be used to distinguish subsets
  of the configuration (e.g., types of particles).
  This selection may be further reduced to single particles.
  These selections are then used to modify a subset of the configuration (e.g.,
  removal and displacement) of a selection of particles/sites.

  The spatial domain contains periodic boundaries and cells.

  Finally, note that an optimization for grand canonical ensemble simulations
  allows for the possibility of ghost particles.
  Ghost particles do not interact and serve as a pool for adding/removing
  particles in an optimized fashion without having to resize the particle
  arrays.
  Thus, one must use caution when accessing particles and sites by indices,
  because these indices may include ghost particles.
 */
class Configuration {
 public:
  //@{
  /** @name Construction
   */

  /**
    args:
    - particle_type[i]: add the i-th type of particle to the configuration.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
      If only one particle type, the "[i]" is optional.
    - wrap: wrap particle centers within domain (default: true).
    - physical_constants: optional class_name of PhysicalConstants.
   */
  explicit Configuration(const argtype& args = argtype());

  /// Same as above, but also set the domain.
  explicit Configuration(std::shared_ptr<Domain> domain,
    const argtype& args = argtype()) : Configuration(args) {domain_ = domain; }

  //@}
  /** @name Typing
    Types of sites and particles.
   */
  //@{

  /// Add a particle type that may exist by LMP file (see FileLMP).
  void add_particle_type(const std::string file_name);

  /// Return the number of particle types.
  int num_particle_types() const { return particle_types_.num(); }

  /// Return the number of site types.
  int num_site_types() const { return unique_types_.num_sites(); }

  /// Return the particle associated with the type.
  const Particle& particle_type(const int type) const {
    return particle_types_.particle(type);
  }

  /// Return the file name used to initialize the particle types.
  std::string type_to_file_name(const int type) const {
    return type_to_file_[type];
  }

  /// Return the particle types.
  const ParticleFactory& particle_types() const { return particle_types_; }

  /// Add a custom type of model parameter.
  /// Name it the same as an atom property before reading file to
  /// make a custom ModelParam.
  void add(std::shared_ptr<ModelParam> param);

  const ModelParams& model_params() const {
    return unique_types_.model_params(); }

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

  /// Add or set model parameter of a given name to value.
  void add_or_set_model_param(const std::string name,
                              const double value) {
    unique_types_.add_or_set_model_param(name, value);
  }

  /// Set the physical constants.
  void set_physical_constants(std::shared_ptr<PhysicalConstants> constants) {
    unique_types_.set_physical_constants(constants); }

  /// Return the physical constants.
  const PhysicalConstants& physical_constants() const {
    return model_params().physical_constants(); }

  /// Return the unique types. Only unique sites and bonds are included.
  /// Thus, the site index is the same as the numeric value for the site type.
  /// And the same for bonds.
  /// This serves as a container for properties based on site or bond type.
  const ParticleFactory& unique_types() const { return unique_types_; }

  /// Return the unique type by individual particle.
  const Particle& unique_type(const int type) const {
    return unique_types_.particle(type);
  }

  /// Return the maximum number of sites in any particle type.
  int max_sites_in_any_particle() const;

  /// Change the site type of a given site in all particles of given type.
  void set_site_type(const int particle_type,
                     const int site,
                     const int site_type);

  //@}
  /** @name Groups
    Groups of sites and particles
   */
  //@{

  /// Add a group (after types are defined but before particles are added).
  void add(std::shared_ptr<Group> group,
    /// Optionally provide a name. If no name is provided, the name is assigned
    /// to be the numerical indices of the order of groups added.
    std::string name = "");

  /// Return the number of group selections.
  int num_groups() const { return static_cast<int>(group_selects_.size()); }

  /// Return the index of the group based on particle types.
  /// If the group does not exist, return -1.
  int particle_type_to_group(const int particle_type) const;

  /// Same as above, except create the group if it does not already exist.
  int particle_type_to_group_create(const int particle_type);

  /// Return the group-based selections.
  const std::vector<Select>& group_selects() const {
    return group_selects_; }

  /// Return the group-based selections by index.
  const Select& group_select(const int index) const {
    return group_selects_[index]; }

  //@}
  /** @name Particles
    Physically existing sites and particles
   */
  //@{

  /// Add a particle of a given type.
  void add_particle_of_type(const int type);

  /// Return particle by index. Note this index is contiguous from values
  /// 0 to num_particles -1, unlike the selection indices (due to ghost
  /// particles).
  /// Note that this method can be slow because the particle index
  /// filters out ghost particles.
  /// Note: this method can be prone to errors if used to define a constant
  /// reference to, for example, site or position in particle.
  const Particle particle(const int index,
    /// Provide a group index to consider only a subset of the configuration.
    /// By default, a value of zero is for the entire configuration.
    const int group = 0) const;

  /// Return selection of all particles and sites in the configuration.
  /// This selection does not include ghost particles.
  const Select& selection_of_all() const { return group_selects_[0]; }

  /// Return the number of particles.
  int num_particles(
    /// Provide a group index to consider only a subset of the configuration.
    /// By default, a value of zero is for the entire configuration.
    const int group = 0) const;

  /// Return the number of sites.
  int num_sites(
    /// Provide a group index as described above.
    const int group = 0) const;

  /// Return the number of particles of a given particle type.
  /// If type == -1, return number of particles of all types.
  int num_particles_of_type(const int type) const;

  /// Return the last particle added.
  const Particle& newest_particle() const {
    return select_particle(newest_particle_index_); }

  /// Return the number of sites of each type in selection.
  std::vector<int> num_sites_of_type(const Select& selection) const;

  /// Return the number of sites of each type in group.
  std::vector<int> num_sites_of_type(const int group_index = 0.) const {
    return num_sites_of_type(group_selects()[group_index]); }

  //@}
  /** @name Modifications
    Modifications to a configuration (e.g., moving, adding or deleting
    particles/sites.
    A subset of the configuration is defined by a Select.
   */
  //@{

  /// Load coordinates by per-site vector containing per-dimension vector.
  /// Requires coordinates for all sites and dimensions.
  void update_positions(const std::vector<std::vector<double> > coords);

  /// Update the positions and properties from a selection.
  void update_positions(const Select& select,
    /// If true, do not wrap. If false, defer to default behavior.
    const bool no_wrap = false,
    /// If true, do not exclude properties.
    const bool no_exclude = false);

  /// Copy the existing particles in a given config by replacing positions.
  /// Types must match config in the same order.
  void copy_particles(const Configuration& config,
    /// Add missing particles of the same type.
    const bool add_missing = false);

  /// Displace selected particle(s). No periodic boundary conditions applied.
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

  /// Set the particle positions of the group to the geometric center of the
  /// site positions.
  void recenter_particle_positions(const int group_index = 0);

  //@{
  /** @name Domain
    A configuration's domain includes periodic boundaries and cells.
   */
  //${

  /// Set the domain.
  void set(std::shared_ptr<Domain> domain);

  /// Return the domain of the configuration.
  const Domain& domain() const { return const_cast<Domain&>(*domain_); }

  // Return the domain.
  Domain * get_domain() { return domain_.get(); }

  /// Set the domain side lengths.
  // HWH consider scaling particles as well
  void set_side_lengths(const Position& sides);

  /// Return the dimensionality of space.
  int dimension() const;

  /// Set whether or not to wrap particles
  void init_wrap(const bool wrap = true) { wrap_ = wrap; }

  //@}
  /** @name Ghosts
    Functions which require knowledge of ghost particles and thus not for
    typical users.
   */
  //@{

  /// Revive the particles in the selection previously removed (ghosts).
  void revive(const Select& selection);

  /// Return the particles.
  /// Warning: typically not for users because it may include ghost particles.
  const ParticleFactory& particles() const { return particles_; }

  /// Return particle by index provided in selection.
  /// Warning: typically not for users because it may include ghost particles.
  const Particle& select_particle(const int index) const {
    return particles_.particle(index); }

  /// Return the selection-based index (includes ghosts) of the last particle
  /// added.
  int newest_particle_index() const { return newest_particle_index_; }

  /// Return ghost particles.
  const std::vector<Select>& ghosts() const { return ghosts_; }

  /// Wrap particle position. The index may include ghost particles.
  void wrap_particle(const int particle_index);

  /// Add a particle of a given type without using ghosts.
  void add_non_ghost_particle_of_type(const int type);

  //@}
  /** @name Sites
    Modify properties of sites directly.
    Note that indices include ghosts.
   */
  //@{

  /// Set selection as physical/nonphysical
  void set_selection_physical(const Select& select, const bool phys);

  /// Set particle property.
  void set_property(const std::string name,
      const double value,
      const int particle_index) {
    particles_.set_property(name, value, particle_index);
  }

  /// Add the property to a site in a particle.
  void add_site_property(const std::string name,
      const double value,
      const int particle_index,
      const int site_index) {
    particles_.add_site_property(name, value, particle_index, site_index);
  }

  /// Add or set the property of a site in a particle.
  void add_or_set_site_property(const std::string name,
      const double value,
      const int particle_index,
      const int site_index) {
    particles_.add_or_set_site_property(name, value,
      particle_index, site_index);
  }

  /// Add or set the property of a site in a particle type.
  void add_or_set_particle_type_site_property(const std::string name,
      const double value,
      const int particle_type,
      const int site_index) {
    particle_types_.add_or_set_site_property(name, value, particle_type,
                                             site_index);
  }

  /// Set the property of a site in a particle by name.
  void set_site_property(const std::string name,
      const double value,
      const int particle_index,
      const int site_index) {
    particles_.set_site_property(name, value, particle_index, site_index);
  }

  /// Set the property of a site in a particle by index.
  void set_site_property(const int index,
      const double value,
      const int particle_index,
      const int site_index) {
    particles_.set_site_property(index, value, particle_index, site_index);
  }

  /// Add an property excluded from Configuration::update_positions()
  void add_excluded_property(const std::string name);

  //@}
  /** @name Checks
    Consistency checks and tests.
   */
  //@{

  /// Check consistency of dimensions and lists.
  void check() const;

  /// Return true if all sites are physical.
  bool are_all_sites_physical() const;

  /// Check if configuration is approximately equivalent.
  /// Not all quantities are checked, including ghosts, etc.
  bool is_equal(const Configuration& configuration,
                const double tolerance) const;

  /// Return the header of the status for periodic output.
  std::string status_header() const;

  /// Return the brief status for periodic output.
  std::string status() const;

  //@}

  // HWH updates entire particle. Optimize by updating per site.
  void synchronize_(const Configuration& config, const Select& perturbed);

  /// Serialize
  void serialize(std::ostream& ostr) const;

  /// Deserialize
  explicit Configuration(std::istream& istr);

 private:
  ParticleFactory particle_types_;
  ParticleFactory unique_types_;
  ParticleFactory particles_;
  std::shared_ptr<Domain> domain_;
  bool wrap_;

  // temporaries (not serialized)
  Arguments args_;
  int newest_particle_index_;
  Select one_site_select_;

  /// Selects based on groups that are continuously updated.
  // HWH currently only updated when adding and removing particles
  // HWH but at some point it should check for positional changes
  // HWH if groups are defined based on positions.
  std::vector<Select> group_selects_;

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
                           const Properties& prop,
                           const std::vector<std::string>& exclude) {
    particles_.replace_properties(particle_index, site_index, prop, exclude); }

  /// Store the excluded properties used in replace_properties_ (optimization).
  std::vector<std::string> excluded_properties_;
  std::vector<std::string> excl_prop_non_usr_ = excluded_properties_;

  /// Update position trackers of a particle (e.g., cell, neighbor, etc).
  void position_tracker_(const int particle_index, const int site_index);

  /// Update position trackers of all sites in a particle.
  void position_tracker_(const int particle_index);

  /// Update position trackers of all particles.
  void position_tracker_();

  /// Add particle to selection.
  void add_to_selection_(const int particle_index,
                         Select * select) const;

  /// Initialize selection based on groups
  void init_selection_(Select * group_select) const;

  /// Remember groups based on types.
  std::vector<int> group_store_particle_type_,
                   group_store_group_index_;

//  /// HWH depreciate one of these.
//  void check_id_(const Select& select) const;
//  void check_id_(const std::string id) const;

  // Ghost particles allow quick addition and deletion of particles for use in
  // the grand canonical ensemble.
  // ghosts are removed from selections and can be brought back by adding.
  // each index represents the particle type.
  std::vector<Select> ghosts_;

  /// Return the number of ghost particles.
  int num_ghosts_() const;

  const Particle& particle_(const int index) {
    return particles_.particle(index);
  }

  /// Store the files used to initialize particle types.
  std::vector<std::string> type_to_file_;

  /// Store the number of particles of each type.
  std::vector<int> num_particles_of_type_;
};

inline std::shared_ptr<Configuration> MakeConfiguration(
    const argtype &args = argtype()
    ) {
  return std::make_shared<Configuration>(args);
}
inline std::shared_ptr<Configuration> MakeConfiguration(
    std::shared_ptr<Domain> domain, const argtype &args = argtype()) {
  return std::make_shared<Configuration>(domain, args);
}

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_CONFIGURATION_H_
