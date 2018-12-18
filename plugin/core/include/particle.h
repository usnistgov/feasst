
#ifndef FEASST_CORE_PARTICLE_H_
#define FEASST_CORE_PARTICLE_H_

#include <vector>
#include "core/include/typed_entity.h"
#include "core/include/site.h"
#include "core/include/domain.h"
#include "core/include/bond.h"

namespace feasst {

/**
  Particles may be a collection of sites (anisotropic) or a single site
  (isotropic).

  The position of a site is often used as its "center" with respect to
  whole-particle operations such as the following:

  1. distance-based particle cut-off methods
  2. wrapping positions about periodic boundary conditions
  3. rigid rotation of the particle

  Each site has its own position which is not relative to the position of
  the particle. Reference sites may be stored optionally.
 */
class Particle : public TypedEntity, public SpatialEntity {
 public:
  /// Add a Site to the Particle.
  void add(Site site) { sites_.push_back(site); }

  /// Return a site for a given index.
  const Site& site(const int index) const { return sites_[index]; }

  /// Set the site of a given index.
  void set_site(const int index, const Site& site) { sites_[index] = site; }

  /// Return the sites.
  const std::vector<Site>& sites() const { return sites_; }

  /// Return the number of Sites.
  int num_sites() const { return sites_.size(); }

  /// Displace the Particle. This also displaces the Sites.
  void displace(const Position& displacement);

  /// Initialize the Particle to be a single site on the origin.
  void default_particle();

  /// Store existing Site(s) for later reference.
  void store_reference_sites() { reference_sites_ = sites_; }

  /// Return true if the particle is isotropic (e.g., 1 or less sites)
  bool is_isotropic();

  /// Check that the dimensionality of the site and particle positions match.
  void check_size();

  /// Increment the types of the sites.
  void increment_site_types(const int increment);

  /// Remove sites and bonds which are of the same type as a previous one.
  void remove_non_unique_types();

  /// Remove a particular site by its index.
  void remove_site(const int index);

  /// Swap the positions of particle with self.
  void replace_position(const Particle& particle);

  /// Replace the position of site by index.
  void replace_position(const int site_index, const Position& replacement);

  /// Update the cells of site index, or the whole particle (default).
  void update_cell(const Cells& cell,
                   const Domain& domain,
                   const int site_index = -1);

  /// Return the number of bonds.
  int num_bonds() const { return static_cast<int>(bonds_.size()); }

  /// Return the bond by index.
  const Bond& bond(const int index) const { return bonds_[index]; }

  /// Return the bonds.
  const std::vector<Bond>& bonds() const { return bonds_; }

  /// Add a bond.
  void add_bond(const Bond bond) { bonds_.push_back(bond); }

  /// Set a bond.
  void set_bond(const int index, const Bond bond) { bonds_[index] = bond; }

  /// Erase all bonds from particle.
  void erase_bonds() { bonds_.clear(); }

//  ~Particle() { check_size(); }
 private:
  std::vector<Site> sites_;
  std::vector<Site> reference_sites_;
  std::vector<Bond> bonds_;

  void update_cell_of_site_(const Cells& cells,
                            const Domain& domain,
                            Site * site);
};

}  // namespace feasst

#endif  // FEASST_CORE_PARTICLE_H_
