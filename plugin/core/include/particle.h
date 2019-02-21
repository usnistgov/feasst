
#ifndef FEASST_CORE_PARTICLE_H_
#define FEASST_CORE_PARTICLE_H_

#include <vector>
#include "core/include/typed_entity.h"
#include "core/include/site.h"
#include "core/include/bond.h"
#include "core/include/debug.h"

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
  /** @name Sites
    Sites of the particle. */
  //@{

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

  /// Remove a particular site by its index.
  void remove_site(const int index);

  /// Store existing Site(s) for later reference.
  void store_reference_sites() { reference_sites_ = sites_; }

  //@}
  /** @name Typing
    Site types and checks
   */
  //@{

  /// Return true if the particle is isotropic (e.g., 1 or less sites)
  bool is_isotropic();

  /// Check that the dimensionality of the site and particle positions match.
  void check_size();

  /// Increment the types of the sites.
  void increment_site_types(const int increment);

  /// Remove sites and bonds which are of the same type as a previous one.
  void remove_non_unique_types();

  //@}
  /** @name Movement
    Move the sites and the particle.
   */
  //@{

  /// Displace the Particle. This also displaces the Sites.
  void displace(const Position& displacement);

  /// Swap the positions of particle with self.
  void replace_position(const Particle& particle);

  /// Replace the position of site by index.
  void replace_position(const int site_index, const Position& replacement);

  /// Return the average of the site positions
  Position average_site_position() const;

  /// Set the particle position as the average of site positions.
  void set_position_as_center();

  //@}
  /** @name Properties
    Change the properties of the sites in the particle
   */
  //@{

  /// Replace the properties of site by index.
  void replace_properties(const int site_index, const Properties& replacement) {
    sites_[site_index].set_properties(replacement); }

  /// Add the property of a site.
  void add_site_property(const std::string name,
      const double value,
      const int site_index) {
    sites_[site_index].add_property(name, value);
  }

  /// Set the property of a site.
  void set_site_property(const std::string name,
      const double value,
      const int site_index) {
    sites_[site_index].set_property(name, value);
  }
  void set_site_property(const int index,
      const double value,
      const int site_index) {
    sites_[site_index].set_property(index, value);
  }

  //@}
  /** @name Bonds
    Bonds between two sites in the particle
   */
  //@{

  /// Return the number of bonds.
  int num_bonds() const { return static_cast<int>(bonds_.size()); }

  /// Return the bond by index.
  const Bond& bond(const int index) const { return bonds_[index]; }

  /// Return the bonds.
  const std::vector<Bond> bonds() const { return bonds_; }

  /// Add a bond.
  void add_bond(const Bond& bond) { bonds_.push_back(bond); }

  /// Set a bond.
  void set_bond(const int index, const Bond& bond) { bonds_[index] = bond; }

  /// Erase all bonds from particle.
  void erase_bonds() { bonds_.clear(); angles_.clear(); }

  //@}
  /** @name Angles
    Angles between three sites in the particle
   */
  //@{

  /// Return the number of angles.
  int num_angles() const { return static_cast<int>(angles_.size()); }

  /// Return the angle by index.
  const Angle& angle(const int index) const { return angles_[index]; }

  /// Return the angles.
  const std::vector<Angle> angles() const { return angles_; }

  /// Add a angle.
  void add_angle(const Angle& angle) { angles_.push_back(angle); }

  /// Set a angle.
  void set_angle(const int index, const Angle& angle) { angles_[index] = angle; }

  //@}
  /** @name Dihedrals
    Dihedrals between four sites in the particle
   */
  //@{

  /// Return the number of dihedrals.
  int num_dihedrals() const { return static_cast<int>(dihedrals_.size()); }

  /// Return the dihedral by index.
  const Dihedral& dihedral(const int index) const { return dihedrals_[index]; }

  /// Return the dihedrals.
  const std::vector<Dihedral> dihedrals() const { return dihedrals_; }

  /// Add a dihedral.
  void add_dihedral(const Dihedral& dihedral) { dihedrals_.push_back(dihedral); }

  /// Set a dihedral.
  void set_dihedral(const int index, const Dihedral& dihedral) { dihedrals_[index] = dihedral; }

  //@}
  /** @name Impropers
    Impropers between four sites in the particle
   */
  //@{

  /// Return the number of impropers.
  int num_impropers() const { return static_cast<int>(impropers_.size()); }

  /// Return the improper by index.
  const Improper& improper(const int index) const { return impropers_[index]; }

  /// Return the impropers.
  const std::vector<Improper> impropers() const { return impropers_; }

  /// Add a improper.
  void add_improper(const Improper& improper) { impropers_.push_back(improper); }

  /// Set a improper.
  void set_improper(const int index, const Improper& improper) { impropers_[index] = improper; }

//  ~Particle() { check_size(); }

 private:
  std::vector<Site> sites_;
  std::vector<Site> reference_sites_;
  std::vector<Bond> bonds_;
  std::vector<Angle> angles_;
  std::vector<Dihedral> dihedrals_;
  std::vector<Improper> impropers_;
};

}  // namespace feasst

#endif  // FEASST_CORE_PARTICLE_H_
