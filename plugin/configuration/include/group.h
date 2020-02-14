
#ifndef FEASST_CONFIGURATION_GROUP_H_
#define FEASST_CONFIGURATION_GROUP_H_

#include <vector>
#include "configuration/include/particle.h"

namespace feasst {

// HWH: chain-setters don't work well with python.
// HWH: use arguments instead
/**
  Define groups based on particle and site types.
  In the future, other metrics may be used, such as position-based ones, etc.
 */
class Group : public PropertiedEntity {
 public:
  Group();

  /// Dynamic groups indicate that they should be updated (default behavior).
  void set_dynamic() { dynamic_ = true; }

  /// Return true if dynamic.
  bool dynamic() const { return dynamic_; }

  /// Add site type as included. Return self for chain setting.
  /// If no site types are defined, then all site types are in the group.
  Group& add_site_type(const int type) {
    site_types_.push_back(type);
    return *this;
  }

  /// Return the list of site types in the group.
  const std::vector<int> site_types() const { return site_types_; }

  /// Add particle type as included. Return self for chain setting.
  /// If no particle types are defined, then all particle types are included.
  Group& add_particle_type(const int type) {
    particle_types_.push_back(type);
    return *this;
  }

  /// Return true if group has no group definitions.
  bool is_empty() const;

  /// Return true if the site is in the group.
  bool is_in(const Site& site) const;

  /// Return true if the particle is in the group.
  bool is_in(const Particle& particle) const;

  /// Remove sites from the particle which are not in the group.
  void remove_sites(Particle * particle) const;

  /// Return the list of site indices in Particle which are in the group.
  std::vector<int> site_indices(const Particle& particle) const;

  /// Return if the group definition is based on location.
  bool is_spatial() const { return spatial_; }

  void serialize(std::ostream& ostr) const;
  Group(std::istream& istr);

 private:
  /// If no site types are listed, do not screen by site types.
  std::vector<int> site_types_;
  std::vector<int> particle_types_;
  bool dynamic_;
  bool spatial_ = false;
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_GROUP_H_
