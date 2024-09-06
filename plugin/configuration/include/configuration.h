
#ifndef FEASST_CONFIGURATION_CONFIGURATION_H_
#define FEASST_CONFIGURATION_CONFIGURATION_H_

#include <memory>
#include <map>
#include <string>
#include <vector>

namespace feasst {

class Domain;
class Group;
class ModelParam;
class ModelParams;
class NeighborCriteria;
class Particle;
class ParticleFactory;
class PhysicalConstants;
class Position;
class Properties;
class Select;
class Site;
class Table3D;
class Table4D;
class Table5D;
class Table6D;

typedef std::map<std::string, std::string> argtype;

/**
  A Configuration contains both the particles and the spatial Domain/boundaries.

  For the particles, this includes both the pool of particles which may
  exist (referred to as particle types) in addition to the particles which
  physically exist (referred to as particles).
  The same is true for site types and sites.

  Groups of different particle/site types and other metrics may be defined.
  These groups then define a selection which can be used to distinguish subsets
  of the configuration (e.g., types of particles).
  This selection may be further reduced to single particles.
  These selections are then used to modify a subset of the configuration (e.g.,
  removal and displacement) of a selection of particles/sites.

  The spatial domain contains periodic boundaries and cells.
 */
class Configuration {
 public:
  //@{
  /** @name Arguments
    The following arguments are parsed in the order listed,
    regardless of the order input by the user.

    - Domain arguments may be parsed here.
    - physical_constants: optional class_name of PhysicalConstants.
      These are typically only used in charged interactions to compute the
      conversion factor between squared charge over distance and energy.
    - particle_type[i]: add the i-th type of particle.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
      If only one particle type, the "[i]" is optional.
    - add_particles_of_type[i]: add this many of i-th type particles.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
    - xyz_file: optionally load FileXYZ if not empty (default: empty).
      Note that Domain tilt factors are not read by FileXYZ.
    - xyz_euler_file: optionally load FileXYZEuler if not empty (default: empty).
    - group[i]: set the name of the "i"-th group.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
      All following arguments of the group are then expected to have the name
      prepended (e.g., "group0 water water_particle_type 0").
    - wrap: wrap particle centers within domain (default: true).
    - [parameter]: optionally set the [parameter] of all types to this value.
      The "[parameter]" is to be substituted for epsilon, sigma, cutoff, etc.
    - [parameter][i]: optionally set the [parameter] of the i-th site type.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
      These are applied after (overriding) the above argument for all types.
    - [parameter][i]_[j]: optionally set the [parameter] of the i-j mixed type.
      The "[i]/[j]" is to be substituted for an integer 0, 1, 2, ...
      These are applied after (overriding) the above argument for single types.
    - [parameter]_mixing_file: override default mixing parameters with this
      file in the three-column space-separated format of "i j [param]_ij".
    - set_cutoff_min_to_sigma: if true and cutoff < sigma, cutoff = sigma
      (default: false). This is typically used for HardSphere models that
      didn't specify cutoff.
   */
  explicit Configuration(argtype args = argtype());
  explicit Configuration(argtype * args);

  //@}
  /** @name Typing
    Types of sites and particles.
   */
  //@{

  /*
    Additional class notes for developers:

    Before particles are physically added to the system, all particle types
    which may exist must be defined.
    Once particles are added to the system, new particle types cannot be added.
    Similarly, site and bond types are also defined at the same time as particle
    types and are referred to as unique types.
    Unique types are stored in essentially the same fashion as particle types,
    except they have been stripped of their non-unique sites and bonds.
    Two different particles cannot share a site type.

    Finally, note that an optimization for grand canonical ensemble simulations
    allows for the possibility of ghost particles.
    Ghost particles do not interact and serve as a pool for adding/removing
    particles in an optimized fashion without having to resize the particle
    arrays.
    Thus, one must use caution when accessing particles and sites by indices,
    because these indices may include ghost particles.
   */

  /// Add a particle type that may exist in the simulation.
  /// See FileParticle.
  void add_particle_type(
    /// Each particle type must have a unique file name.
    const std::string file_name,
    /// Optionally, append to name to use same file but keep unique names.
    const std::string append = "");

  /// Return the number of particle types.
  int num_particle_types() const;

  /// Return the number of site types.
  int num_site_types() const;

  /// Return the number of bond types.
  int num_bond_types() const;

  /// Return the number of angle types.
  int num_angle_types() const;

  /// Return the number of dihedral types.
  int num_dihedral_types() const;

  /// Return the particle associated with the type.
  const Particle& particle_type(const int type) const;

  /// Return the file name used to initialize the particle types.
  const std::string& type_to_file_name(const int type) const {
    return type_to_file_[type]; }

  /// Return the particle types.
  const ParticleFactory& particle_types() const;

  /// Add a custom type of model parameter.
  /// Name it the same as an atom property before reading file to
  /// make a custom ModelParam.
  void add(std::shared_ptr<ModelParam> param);

  /// Return the model parameters (e.g., sigma, epsilon, etc).
  const ModelParams& model_params() const;

  /// Modify model parameter of a given site type and name to value.
  void set_model_param(const std::string name,
                       const int site_type,
                       const double value);

  /// Modify a mixed model parameter of given site types and name to value.
  void set_model_param(const std::string name,
                       const int site_type1,
                       const int site_type2,
                       const double value);

  /// Set mixed model parameters using a file.
  void set_model_param(const std::string name, const std::string filename);

  /// Add model parameter of a given name to value.
  void add_model_param(const std::string name, const double value);

  /// Add or set model parameter of a given name to value.
  void add_or_set_model_param(const std::string name, const double value);

  /// Set the physical constants.
  void set_physical_constants(std::shared_ptr<PhysicalConstants> constants);

  /// Return the physical constants.
  const PhysicalConstants& physical_constants() const;

  /// Return the unique types. Only unique sites and bonds are included.
  /// Thus, the site index is the same as the numeric value for the site type.
  /// And the same for bonds.
  /// This serves as a container for properties based on site or bond type.
  const ParticleFactory& unique_types() const;

  /// Return the unique type by individual particle.
  const Particle& unique_type(const int type) const;

  /// Return the site of unique type by individual particle and site type.
  const Site& unique_type(const int ptype, const int stype) const;

  /// Return the maximum number of sites in any particle type.
  int max_sites_in_any_particle() const;

  /// Change the site type of a given site in all particles of given type.
  void set_site_type(const int particle_type,
                     const int site,
                     const int site_type);

  /// Return the particle type of a given site type.
  int site_type_to_particle_type(const int site_type) const;

  /// Return, for each particle type, the number of sites of each type.
  std::vector<std::vector<int> > num_site_types_per_particle_type() const;

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
  int num_groups() const;

  /// Return the index of the group based on particle types.
  /// If the group does not exist, return -1.
  int particle_type_to_group(const int particle_type) const;

  /// Same as above, except create the group if it does not already exist.
  int particle_type_to_group_create(const int particle_type);

  /// Return the index of the group with the given name.
  int group_index(const std::string& name) const;

  /// Return the group-based selections.
  const std::vector<std::shared_ptr<Select> >& group_selects() const;

  /// Return the group-based selections by index.
  const Select& group_select(const int index) const;

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
  const Particle& particle(const int index) const;

  /// Same as above, but filter the particle by the group index.
  /// Returns a copy of the filtered particle, instead of a constant reference.
  /// Note: this method can be prone to errors if used to define a constant
  /// reference to, for example, site or position in particle.
  Particle particle(const int index,
    /// Provide a group index to consider only a subset of the configuration.
    const int group) const;

  /// Return selection of all particles and sites in the configuration.
  /// This selection does not include ghost particles.
  const Select& selection_of_all() const;

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

  /// Same as above, but optimized to use existing data structure.
  void num_sites_of_type(const Select& selection, std::vector<int> * num) const;

  /// Return the number of sites of each type in group.
  std::vector<int> num_sites_of_type(const int group_index = 0.) const {
    return num_sites_of_type(*group_selects()[group_index]); }

  /// Same as above, but optimized to use existing data structure.
  void num_sites_of_type(const int group_index, std::vector<int> * num) const {
    num_sites_of_type(*group_selects()[group_index], num); }

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

  /**
    Load coordinates and orientations with a per-site vector containing
    per-dimension vector.
    Requires coordinates and orientations for all sites and dimensions.
   */
  void update_positions(const std::vector<std::vector<double> > coords,
                        const std::vector<std::vector<double> > eulers);

  /// Update the positions and properties from a selection.
  /// Includes euler angles.
  void update_positions(const Select& select,
    /// If true, do not wrap. If false, defer to default behavior.
    const bool no_wrap = false);

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

  /// Return the number of cell lists.
  int num_cell_lists() const { return num_cell_lists_; }

  /// Increment the number of cell lists.
  void increment_num_cell_lists() { ++num_cell_lists_; }

  /*
    Change the volume. Also, update cells.

    args:
    - dimension: index of dimension to change. If -1, change all (default: -1)
    - scale_particles: if true, scale particle positions (default: true).
   */
  void change_volume(const double delta_volume, argtype args);
  void change_volume(const double delta_volume, argtype * args);

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
  const ParticleFactory& particles() const;

  // Return the particles.
  ParticleFactory * get_particles_();

  /// Return particle by index provided in selection.
  /// Warning: typically not for users because it may include ghost particles.
  const Particle& select_particle(const int index) const;

  /// Return the selection-based index (includes ghosts) of the last particle
  /// added.
  int newest_particle_index() const { return newest_particle_index_; }

  /// Return ghost particles.
  const std::vector<std::shared_ptr<Select> >& ghosts() const;

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
  void set_property(const std::string name, const double value,
    const int particle_index);

  /// Add the property to a site in a particle.
  void add_site_property(const std::string name, const double value,
    const int particle_index,  const int site_index);

  /// Add or set the property of a site in a particle.
  void add_or_set_site_property(const std::string name, const double value,
      const int particle_index, const int site_index);

  /// Add or set the property of a site in a particle type.
  void add_or_set_particle_type_site_property(const std::string name,
      const double value,
      const int particle_type,
      const int site_index);

  /// Set the property of a site in a particle by name.
  void set_site_property(const std::string name, const double value,
    const int particle_index, const int site_index);

  /// Set the property of a site in a particle by index.
  void set_site_property(const int index, const double value,
    const int particle_index, const int site_index);

  /// Change select to a given particle type.
  void set_particle_type(const int particle_type,
                         const Select& select);

  //@}
  /** @name Neighbor Criteria
    Define and store various criteria used for defining neighbors
   */
  //@{

  /// Add a NeighborCriteria.
  void add(std::shared_ptr<NeighborCriteria> neighbor_criteria);

  /// Return a NeighborCriteria by index in order added.
  const NeighborCriteria& neighbor_criteria(const int index) const;

  /// Return a NeighborCriteria by index in order added.
  const std::vector<std::shared_ptr<NeighborCriteria> >& neighbor_criteria()
    const;

  // Return a NeighborCriteria by index in order added.
  NeighborCriteria * get_neighbor_criteria(const int index);

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

  /// Check that the dimensions of all particle types and domain match
  void check_dimensions() const;

  /// Return the header of the status for periodic output.
  std::string status_header(const std::string append = "") const;

  /// Return the brief status for periodic output.
  std::string status() const;

  /// Return the configuration in human readable format for testing.
  std::string str() const;

  //@}
  /** @name Tables
    Shared tables that can be used in multiple VisitModel.
   */
  //@{

  const std::vector<std::vector<std::shared_ptr<Table3D> > >& table3d() const;
  const std::vector<std::vector<std::shared_ptr<Table4D> > >& table4d() const;
  const std::vector<std::vector<std::shared_ptr<Table5D> > >& table5d() const;
  const std::vector<std::vector<std::shared_ptr<Table6D> > >& table6d() const;
  std::vector<std::vector<std::shared_ptr<Table3D> > > * get_table3d();
  std::vector<std::vector<std::shared_ptr<Table4D> > > * get_table4d();
  std::vector<std::vector<std::shared_ptr<Table5D> > > * get_table5d();
  std::vector<std::vector<std::shared_ptr<Table6D> > > * get_table6d();

  //@}

  // HWH updates entire particle. Optimize by updating per site.
  void synchronize_(const Configuration& config, const Select& perturbed);

  /// Serialize
  void serialize(std::ostream& ostr) const;

  /// Deserialize
  explicit Configuration(std::istream& istr);

 private:
  std::shared_ptr<ParticleFactory> particle_types_;
  std::shared_ptr<ParticleFactory> unique_types_;
  std::shared_ptr<ParticleFactory> particles_;
  std::shared_ptr<Domain> domain_;
  bool wrap_;
  int num_cell_lists_ = 0;

  // temporaries (not serialized)
  int newest_particle_index_;

  /// Selects based on groups that are continuously updated.
  // HWH currently only updated when adding and removing particles
  // HWH but at some point it should check for positional changes
  // HWH if groups are defined based on positions.
  std::vector<std::shared_ptr<Select> > group_selects_;

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

  /// Replace properties of site in particle.
  void replace_properties_(const int particle_index,
                           const int site_index,
                           const Properties& prop);

  /// Update position trackers of a particle (e.g., cell, neighbor, etc).
  void position_tracker_(const int particle_index, const int site_index);

  /// Update position trackers of all sites in a particle.
  void position_tracker_(const int particle_index);

  /// Update position trackers of all particles.
  void position_tracker_();

  /// Add particle to selection.
  void add_to_selection_(const int particle_index,
                         Select * select) const;

  /// Update particle in selection.
  void update_selection_(const int particle_index,
                         Select * select) const;

  /// Initialize selection based on groups
  /// Initialize selection based on groups
  void init_selection_(Select * group_select) const;

  /// Remember groups based on types.
  std::vector<int> group_store_particle_type_,
                   group_store_group_index_;

//  /// HWH deprecate one of these.
//  void check_id_(const Select& select) const;
//  void check_id_(const std::string id) const;

  // Ghost particles allow quick addition and deletion of particles for use in
  // the grand canonical ensemble.
  // ghosts are removed from selections and can be brought back by adding.
  // each index represents the particle type.
  std::vector<std::shared_ptr<Select> > ghosts_;

  /// Return the number of ghost particles.
  int num_ghosts_() const;

  const Particle& particle_(const int index);

  /// Store the files used to initialize particle types.
  std::vector<std::string> type_to_file_;

  /// Store the number of particles of each type.
  std::vector<int> num_particles_of_type_;
  std::vector<std::shared_ptr<NeighborCriteria> > neighbor_criteria_;

  // Do not serialize possibly huge tables.
  std::vector<std::vector<std::shared_ptr<Table3D> > > table3d_;
  std::vector<std::vector<std::shared_ptr<Table4D> > > table4d_;
  std::vector<std::vector<std::shared_ptr<Table5D> > > table5d_;
  std::vector<std::vector<std::shared_ptr<Table6D> > > table6d_;
};

inline std::shared_ptr<Configuration> MakeConfiguration(
    argtype args = argtype()) {
  return std::make_shared<Configuration>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_CONFIGURATION_H_
