
#ifndef FEASST_CONFIGURATION_PARTICLE_FACTORY_H_
#define FEASST_CONFIGURATION_PARTICLE_FACTORY_H_

#include <memory>
#include <string>
#include <vector>
#include "configuration/include/particle.h"
#include "configuration/include/model_params.h"
#include "configuration/include/group.h"

namespace feasst {

/**
  A container for a list of particles.
  There are effectively three different uses:

    1. A list of particles that exist in a simulation. In this case there can
       be multiple particles and sites of the same type.

    2. A list of particle types that may exist in a simulation. In this case
       there can not be multiple particles of the same type.
       This is enforced by ParticleFactory::unique_particles().

    3. A list of site and bond types (contained within particles) that may exist
       in a simulation.
       In this case, there can not be multiple sites or bonds of the same type.
       This contains site-based and bond-based properties.
       This is enforced by ParticleFactory::unique_types().
 */
class ParticleFactory {
 public:
  ParticleFactory() {}

  /// Adjust the site types of added particles to ensure uniqueness.
  /// Returns self for chain setting.
  ParticleFactory& unique_particles();

  /**
    Only add sites or bonds which are new as a holder for site- or bond- type
    properties.
    Note that ParticleFactory::unique_particles() is also applied.
    Returns self for chain setting.
   */
  ParticleFactory& unique_types();

  /**
    Check that no site, particle, bond or angle type is skipped.
    As particles and sites are added, a never before seen type must be only
    one higher than the previous maximum index.
    Types begin with 0.

    For example, a particle with sites of type {0, 2} is invalid because
    the second site of type "2" appears without a previous type "1" site.

    In another example, a particle with sites of type {1,2} is invalid because
    the site type doesn't begin with 0.

    Note: for partial configurations by groups, these checks may not apply.
   */
  void check_types(int * num_site_types, int * num_particle_types,
                   int * num_bond_types, int * num_angle_types) const;

  /// Check that no site, particle, bond or angle type is skipped.
  void check_types() const;

  /// Check that no site type is skipped. Return number of site types.
  int check_site_types() const;

  /// Check that no particle type is skipped. Return number of particle types.
  int check_particle_types() const;

  /// Check that no bond type is skipped. Return number of bond types.
  int check_bond_types() const;

  /// Check that no angle type is skipped. Return number of angle types.
  int check_angle_types() const;

  /// Add particle by file.
  void add(const std::string file_name);

  /// Add a particle.
  void add(const Particle& particle);

  /// Remove particle by index.
  void remove(const int particle_index);

  /// Return particle by index.
  const Particle& particle(const int particle_index) const {
    return particles_[particle_index]; }

  /// Return particles.
  const std::vector<Particle>& particles() const { return particles_; }

  /// Return the number of particles.
  int num() const { return static_cast<int>(particles_.size()); }

  /// Replace position of the particle by index.
  void replace_position(const int particle_index, const Particle& replacement);

  /// Replace position of the site by index.
  void replace_position(const int particle_index,
                        const int site_index,
                        const Position& replacement);

  /// Replace properties of the site by index.
  void replace_properties(const int particle_index,
                          const int site_index,
                          const Properties& replacement) {
    particles_[particle_index].replace_properties(site_index, replacement); }

  /// Scale particle positions by a constant factor in the given dimension.
  void scale_particle_positions(
    /// Scale in this dimension. If -1, scale in all dimensions.
    const int dimen,
    const double factor);

  /// Check consistency of dimensionality of positions of particles and sites.
  /// By default, for dimension == -1, determine automatically.
  void check(const int dimension = -1) const;

  /// Return the number of site types.
  int num_site_types() const { return check_site_types(); }

  /// Return the number of sites.
  int num_sites() const;

  /// Return the number of particle types.
  int num_particle_types() const { return check_particle_types(); }

  /// Return the number of bond types.
  int num_bond_types() const { return check_bond_types(); }

  /// Return the number of bonds.
  int num_bonds() const;

  /// Return the number of angle types.
  int num_angle_types() const { return check_angle_types(); }

  /// Return the number of angles.
  int num_angles() const;

  /// Change the site type of a given site in a particle.
  void set_site_type(const int particle,
                     const int site,
                     const int site_type) {
    particles_[particle].set_site_type(site, site_type); }

  /// Remove particles and sites based on the group.
  void remove(const Group group);

  /// Displace the particle with given index.
  void displace(const int particle_index, const Position& displacement) {
    particles_[particle_index].displace(displacement); }

  /// Return the model parameters.
  const ModelParams& model_params() const { return model_params_; }

  /// Add a custom type of model parameter.
  void add(const std::shared_ptr<ModelParam> param) {
    model_params_.add(param); }

  /// Modify model parameter of a given site type and name to value.
  void set_model_param(const std::string name,
                       const int site_type,
                       const double value) {
    model_params_.set(name, site_type, value); }

  /// Modify a mixed model parameter of given site types and name to value.
  void set_model_param(const char* name,
                       const int site_type1,
                       const int site_type2,
                       const double value) {
    model_params_.set(name, site_type1, site_type2, value); }

  /// Add model parameter of a given name to value.
  void add_model_param(const std::string name,
                       const double value) {
    model_params_.add_property(name, value); }

  /// Add or set model parameter of a given name to value.
  void add_or_set_model_param(const std::string name,
                              const double value) {
    model_params_.add_or_set_property(name, value); }

  /// Set the minimum cutoff to sigma.
  void set_cutoff_min_to_sigma() { model_params_.set_cutoff_min_to_sigma(); }

  /// Set the physical constants in model parameters.
  void set_physical_constants(std::shared_ptr<PhysicalConstants> constants) {
    model_params_.set_physical_constants(constants); }

  /// Set site as physical/nonphysical.
  void set_site_physical(const int particle, const int site, const bool phys) {
    particles_[particle].set_site_physical(site, phys); }

  /// Add particle property.
  void add_property(const std::string name,
      const double value,
      const int particle_index) {
    particles_[particle_index].add_property(name, value); }

  /// Add property to all particles.
  void add_property(const std::string name, const double value) {
    for (Particle& part : particles_) {
      part.add_property(name, value);
    }
  }

  /// Set particle property.
  void set_property(const std::string name,
      const double value,
      const int particle_index) {
    particles_[particle_index].set_property(name, value); }

  /// Add the property of sites in a particle.
  void add_site_property(const std::string name,
      const double value,
      const int particle_index,
      const int site_index) {
    particles_[particle_index].add_site_property(name, value, site_index); }

  /// Add or set the property of sites in a particle.
  void add_or_set_site_property(const std::string name,
      const double value,
      const int particle_index,
      const int site_index) {
    particles_[particle_index].add_or_set_site_property(name, value,
                                                        site_index); }

  /// Set the property of sites in a particle.
  void set_site_property(const std::string name,
      const double value,
      const int particle_index,
      const int site_index) {
    particles_[particle_index].set_site_property(name, value, site_index); }

  /// Set the property of sites in a particle by index instead of name.
  void set_site_property(const int index,
      const double value,
      const int particle_index,
      const int site_index) {
    particles_[particle_index].set_site_property(index, value, site_index); }

  // This interface is for optimization and not for typical use
  Particle * get_particle(const int index) { return &particles_[index]; }

  void serialize(std::ostream& ostr) const;
  explicit ParticleFactory(std::istream& istr);

 private:
  std::vector<Particle> particles_;
  bool unique_particles_ = false;
  bool unique_types_ = false;
  ModelParams model_params_;
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_PARTICLE_FACTORY_H_
