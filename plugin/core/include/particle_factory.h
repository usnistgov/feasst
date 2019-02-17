
#ifndef FEASST_CORE_PARTICLE_FACTORY_H_
#define FEASST_CORE_PARTICLE_FACTORY_H_

#include <vector>
#include "core/include/particle.h"
#include "core/include/model_params.h"
#include "core/include/group.h"

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
       In this case, there can not be multiple sites or bonds of the same
       type. This contains site-based and bond-based properties.
       This is enforced by ParticleFactory::unique_types().
 */
class ParticleFactory {
 public:
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
    Check that no site or particle type is skipped.
    As particles and sites are added, a never before seen type must be only
    one higher than the previous maximum index.
    Types begin with 0.

    For example, a particle with sites of type {0, 2} is invalid because
    the second site of type "2" appears without a previous type "1" site.

    In another example, a particle with sites of type {1,2} is invalid because
    the site type doesn't begin with 0.

    Note: for partial configurations by groups, these checks may not apply.
   */
  void check_types(int * num_site_types, int * num_particle_types) const;

  /// Check that no site or particle type is skipped.
  void check_types() const;

  /// Check that no site type is skipped.
  /// \return number of site types.
  int check_site_types() const;

  /// Check that no particle type is skipped.
  /// \return number of particle types.
  int check_particle_types() const;

  /// Add particle by file.
  void add(const char* file_name);

  /// Add a particle.
  void add(const Particle& particle);

  /// Remove particle by index.
  void remove(const int particle_index);

  /// Return particle by index.
  const Particle& particle(const int particle_index) const;

  /// Return particles.
  std::vector<Particle> particles() const { return particles_; }

  /// Return the number of particles.
  int num() const { return static_cast<int>(particles_.size()); }

  /// Replace position of the particle by index.
  void replace_position(const int particle_index, const Particle& replacement);

  /// Replace position of the site by index.
  void replace_position(const int particle_index,
                        const int site_index,
                        const Position& replacement);

  /// Replace position of particle but not site.
  void replace_position(const int particle_index,
                        const Position& replacement);

  /// Replace properties of the site by index.
  void replace_properties(const int particle_index,
                          const int site_index,
                          const Properties& replacement) {
    particles_[particle_index].replace_properties(site_index, replacement); }

  /// Check consistency of dimensionality of positions of particles and sites.
  /// By default, for dimension == -1, determine automatically.
  void check_size(const int dimension = -1) const;

  /// Return the number of site types.
  int num_site_types() const { return check_site_types(); }

  /// Return the number of sites.
  int num_sites() const;

  /// Return the number of particle types.
  int num_particle_types() const { return check_particle_types(); }

  /// Remove particles and sites based on the group.
  void remove(const Group group);

  /// Displace the particle with given index.
  void displace(const int particle_index, const Position& displacement) {
    particles_[particle_index].displace(displacement);
  }

  /// Return the model parameters.
  const ModelParams& model_params() const { return model_params_; }

  /// Modify model parameter of a given site type and name to value.
  void set_model_param(const char* name,
                       const int site_type,
                       const double value) {
    model_params_.set(name, site_type, value);
  }

  /// Add model parameter of a given name to value.
  void add_model_param(const std::string name,
                       const double value) {
    model_params_.add_property(name, value);
  }

  /// Add the property of sites in a particle.
  void add_site_property(const std::string name,
      const double value,
      const int particle_index,
      const int site_index) {
    particles_[particle_index].add_site_property(name, value, site_index);
  }

  /// Set the property of sites in a particle.
  void set_site_property(const std::string name,
      const double value,
      const int particle_index,
      const int site_index) {
    particles_[particle_index].set_site_property(name, value, site_index);
  }
  void set_site_property(const int index,
      const double value,
      const int particle_index,
      const int site_index) {
    particles_[particle_index].set_site_property(index, value, site_index);
  }

 private:
  std::vector<Particle> particles_;
  bool unique_particles_ = false;
  bool unique_types_ = false;
  ModelParams model_params_;
};

}  // namespace feasst

#endif  // FEASST_CORE_PARTICLE_FACTORY_H_
