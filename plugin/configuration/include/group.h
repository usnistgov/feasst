
#ifndef FEASST_CONFIGURATION_GROUP_H_
#define FEASST_CONFIGURATION_GROUP_H_

#include <memory>
#include <vector>
#include <map>
#include <string>
#include "configuration/include/properties.h"

namespace feasst {

class Particle;
class ParticleFactory;
class Site;

typedef std::map<std::string, std::string> argtype;

// HWH In the future, other metrics may be used, such as position-based, etc.
/**
  Define groups based on particle and site types.
 */
class Group : public PropertiedEntity {
 public:
  //@{
  /** @name Arguments
    - prepend: expect all other arguments to have this prepended with underscore.
      For example, if prepend==water, the following argument would expect
      "water_site_type0" (default: empty).
    - site_type: add site type(s). If none, all sites included.
      Multiple can be provided as comma-separated values.
    - particle_type: add particle type(s). If none, all included.
      Multiple can be provided as comma-separated values.
    - particle_index: add particle index(es). If none, all included.
      Multiple can be provided as comma-separated values.
    - dynamic: set true if groups should be updated (default: true).
    - spatial: set true if group is based on location (default: false).
   */
  explicit Group(argtype args = argtype());
  explicit Group(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Convert site type names to interger indexes for optimization.
  /// Otherwise, assumes name is index (and index is -1 if stoi fails).
  void name_to_index(const ParticleFactory& unique_types);

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

  //@}
 private:
  /// If no types or indices are listed, do not screen by types or indices.
  std::vector<std::string> site_type_names_;
  std::vector<int> site_types_;
  std::vector<std::string> particle_type_names_;
  std::vector<int> particle_types_;
  // std::vector<int> site_indices_;
  std::vector<int> particle_indices_;
  bool dynamic_;
  bool spatial_;
};

inline std::shared_ptr<Group> MakeGroup(argtype args = argtype()) {
  return std::make_shared<Group>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_GROUP_H_
