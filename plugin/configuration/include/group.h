
#ifndef FEASST_CONFIGURATION_GROUP_H_
#define FEASST_CONFIGURATION_GROUP_H_

#include <memory>
#include <vector>
#include "utils/include/arguments.h"
#include "configuration/include/particle.h"

namespace feasst {

/**
  Define groups based on particle and site types.
  In the future, other metrics may be used, such as position-based ones, etc.
 */
class Group : public PropertiedEntity {
 public:
  /**
    args:
    - prepend: expect all other arguments to have this prepended with underscore.
      For example, if prepend==water, the following argument would expect
      "water_site_type0" (default: empty).
    - site_type[i]: add the i-th site type. If none, all sites included.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
      If only one site type, the "[i]" is optional.
    - particle_type[i]: add the i-th particle type. If none, all included.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
      If only one particle type, the "[i]" is optional.
    - particle_index[i]: add the i-th particle index. If none, all included.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
      If only one particle index, the "[i]" is optional.
    - dynamic: set true if groups should be updated (default: true).
    - spatial: set true if group is based on location (default: false).
   */
  explicit Group(argtype args = argtype());
  explicit Group(argtype * args);

  /// Return the list of site types in the group.
  const std::vector<int> site_types() const { return site_types_; }

  /// Return the list of particle types in the group.
  const std::vector<int> particle_types() const { return particle_types_; }

  /// Return true if dynamic.
  bool is_dynamic() const { return dynamic_; }

  /// Return if the group definition is based on location.
  bool is_spatial() const { return spatial_; }

  /// Return true if group has no group definitions.
  bool is_empty() const;

  /// Return true if the site is in the group.
  bool is_in(const Site& site) const;

  /// Return true if the particle is in the group.
  bool is_in(const Particle& particle, const int particle_index) const;

  /// Remove sites from the particle which are not in the group.
  void remove_sites(Particle * particle) const;

  /// Return the list of site indices in Particle which are in the group.
  std::vector<int> site_indices(const Particle& particle) const;

  void serialize(std::ostream& ostr) const;
  explicit Group(std::istream& istr);

 private:
  /// If no types or indices are listed, do not screen by types or indices.
  std::vector<int> site_types_;
  std::vector<int> particle_types_;
  //std::vector<int> site_indices_;
  std::vector<int> particle_indices_;
  bool dynamic_;
  bool spatial_;
};

inline std::shared_ptr<Group> MakeGroup(argtype args = argtype()) {
  return std::make_shared<Group>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_GROUP_H_
