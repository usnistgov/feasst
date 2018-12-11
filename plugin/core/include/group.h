
#ifndef FEASST_CORE_GROUP_H_
#define FEASST_CORE_GROUP_H_

#include <vector>
#include "core/include/particle.h"

namespace feasst {

/**
  Define groups based on particle and site types.
  In the future, other metrics may be used, such as position-based ones, etc.
 */
class Group {
 public:
  Group();

  /// Dynamic groups indicate that they should be updated (default behavior).
  void set_dynamic() { dynamic_ = true; }
  bool dynamic() const { return dynamic_; }

  /// Add site type as included. Return self for chain setting.
  /// If no site types are defined, then all site types are in the group.
  Group& add_site_type(const int type) {
    site_types_.push_back(type);
    return *this;
  }

  /// Add particle type as included. Return self for chain setting.
  /// If no particle types are defined, then all particle types are included.
  Group& add_particle_type(const int type) {
    particle_types_.push_back(type);
    return *this;
  }

  /// Return true if group has no group definitions.
  bool empty() const;

  /// Return true if the site is in the group.
  bool is_in(const Site& site) const;

  /// Return true if the particle is in the group.
  bool is_in(const Particle& particle) const;

  /// Remove sites from the particle which are not in the group.
  void remove_sites(Particle * particle,
                    std::vector<int> * full_to_partial_site = NULL,
                    std::vector<int> * partial_to_full_site = NULL) const;

  /// Return particle with sites removed as described above.
  Particle remove_sites(const Particle particle,
                        std::vector<int> * full_to_partial_site = NULL,
                        std::vector<int> * partial_to_full_site = NULL) const;

 private:
  /// If no site types are listed, do not screen by site types.
  std::vector<int> site_types_;
  std::vector<int> particle_types_;
  bool dynamic_;
};

}  // namespace feasst

#endif  // FEASST_CORE_GROUP_H_
