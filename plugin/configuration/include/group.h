
#ifndef FEASST_CONFIGURATION_GROUP_H_
#define FEASST_CONFIGURATION_GROUP_H_

#include <vector>
#include "utils/include/arguments.h"
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
  /**
    args:
    - add_site_type: add a site type. If none, all sites included.
    - add_particle_type: add a particle type. If none, all particles included.
    - dynamic: set true if groups should be updated (default: true).
    - spatial: set true if group is based on location (default: false).
   */
  explicit Group(const argtype& args = argtype());

  /// Return the list of site types in the group.
  const std::vector<int> site_types() const { return site_types_; }

  /// Return true if dynamic.
  bool is_dynamic() const { return dynamic_; }

  /// Return if the group definition is based on location.
  bool is_spatial() const { return spatial_; }

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

  void serialize(std::ostream& ostr) const;
  Group(std::istream& istr);

 private:
  /// If no site types are listed, do not screen by site types.
  std::vector<int> site_types_;
  std::vector<int> particle_types_;
  bool dynamic_;
  bool spatial_;
};

inline std::shared_ptr<Group> MakeGroup(const argtype &args = argtype()) {
  return std::make_shared<Group>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_GROUP_H_
