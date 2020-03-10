
#ifndef FEASST_CONFIGURATION_PARTICLE_H_
#define FEASST_CONFIGURATION_PARTICLE_H_

#include <vector>
#include <string>
#include "configuration/include/properties.h"
#include "configuration/include/typed_entity.h"
#include "configuration/include/site.h"
#include "configuration/include/bond.h"
#include "utils/include/debug.h"

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
class Particle : public PropertiedEntity,
                 public TypedEntity,
                 public SpatialEntity {
 public:
  Particle() : PropertiedEntity(), TypedEntity(), SpatialEntity() {}

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

  //@}
  /** @name Typing
    Site types and checks
   */
  //@{

  /// Return true if the particle is isotropic (e.g., 1 or less sites)
  bool is_isotropic();

  /// Set a site as physical/nonphysical.
  void set_site_physical(const int site, const bool physical) {
    sites_[site].set_physical(physical); }

  /// Check that the dimensionality of the site and particle positions match.
  void check();

  /// Increment the types of the sites.
  void increment_site_types(const int increment);

  /// Remove sites and bonds which are of the same type as a previous one.
  void remove_non_unique_types();

  /// Set the site type.
  void set_site_type(const int site, const int type) {
    sites_[site].set_type(type); }

  int num_sites_of_type(const int type) const;

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
  void replace_properties(const int site_index,
      const Properties& replacement,
      const std::vector<std::string>& exclude) {
    sites_[site_index].set_properties(replacement, exclude); }

  /// Add the property of a site.
  void add_site_property(const std::string name,
      const double value,
      const int site_index) {
    sites_[site_index].add_property(name, value);
  }

  /// Add or set the property of a site.
  void add_or_set_site_property(const std::string name,
      const double value,
      const int site_index) {
    sites_[site_index].add_or_set_property(name, value);
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
  void add_bond(const Bond& bond);

  /// Add a property to a bond.
  void add_bond_property(const int bond, const std::string name,
    const double value) { bonds_[bond].add_property(name, value); }

  /// Erase all bonds from particle.
  void erase_bonds();

  /// Find the bond between the given site indices.
  const Bond& bond(const int site_index1, const int site_index2) const;

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

  /// Add an angle.
  void add_angle(const Angle& angle);

  /// Add a property to a bond.
  void add_angle_property(const int angle, const std::string name,
    const double value) { angles_[angle].add_property(name, value); }

  /// Find the angle between given sites 1-2-3, with 2 as the vertex.
  const Angle& angle(const int site_index1,
                     const int site_index2,
                     const int site_index3) const;

//  //@}
//  /** @name Dihedrals
//    Dihedrals between four sites in the particle
//   */
//  //@{
//
//  /// Return the number of dihedrals.
//  int num_dihedrals() const { return static_cast<int>(dihedrals_.size()); }
//
//  /// Return the dihedral by index.
//  const Dihedral& dihedral(const int index) const {
//    return dihedrals_[index]; }
//
//  /// Return the dihedrals.
//  const std::vector<Dihedral> dihedrals() const { return dihedrals_; }
//
//  /// Add a dihedral.
//  void add_dihedral(const Dihedral& dihedral) {
//    dihedrals_.push_back(dihedral); }
//
//  /// Set a dihedral.
//  void set_dihedral(const int index, const Dihedral& dihedral) {
//    dihedrals_[index] = dihedral; }
//
//  //@}
//  /** @name Impropers
//    Impropers between four sites in the particle
//   */
//  //@{
//
//  /// Return the number of impropers.
//  int num_impropers() const { return static_cast<int>(impropers_.size()); }
//
//  /// Return the improper by index.
//  const Improper& improper(const int index) const {
//    return impropers_[index]; }
//
//  /// Return the impropers.
//  const std::vector<Improper> impropers() const { return impropers_; }
//
//  /// Add an improper.
//  void add_improper(const Improper& improper) {
//    impropers_.push_back(improper); }
//
//  /// Set an improper.
//  void set_improper(const int index, const Improper& improper) {
//    impropers_[index] = improper; }

//  ~Particle() { check(); }

  const std::vector<std::vector<int> >& bond_list() const { return bond_list_; }
  const std::vector<std::vector<int> >& bond_neighbor() const {
    return bond_neighbor_; }
  const std::vector<std::vector<int> >& angle_list() const {
    return angle_list_; }

  void serialize(std::ostream& ostr) const;
  explicit Particle(std::istream& istr);

 private:
  std::vector<Site> sites_;
  std::vector<Bond> bonds_;
  std::vector<Angle> angles_;
//  std::vector<Dihedral> dihedrals_;
//  std::vector<Improper> impropers_;

  // the first dimension is the site index, the second is:
  // bond index
  std::vector<std::vector<int> > bond_list_;
  // the index of the other site
  std::vector<std::vector<int> > bond_neighbor_;
  // angle index
  std::vector<std::vector<int> > angle_list_;

  void resize_list_(std::vector<std::vector<int> > * list);
  void add_bond_(const Bond& bond, const int index,
    std::vector<std::vector<int> > * list);
  void add_bond_neighbor_(const Bond& bond,
                          std::vector<std::vector<int> > * list);
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_PARTICLE_H_
