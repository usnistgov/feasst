
#ifndef FEASST_CONFIGURATION_BOND_H_
#define FEASST_CONFIGURATION_BOND_H_

#include <string>
#include <vector>
#include "configuration/include/typed_entity.h"
#include "configuration/include/properties.h"

namespace feasst {

/**
  Sites within the same particle may be bonded.
  The indices of the sites which are bonded are stored here.
  The type of the bond is used to determine the bond model.
 */
class Bond : public PropertiedEntity, public TypedEntity {
 public:
  Bond() { set_name("bond"); }

  /// Return the indices of the sites involved in the bond within a particle.
  std::vector<int> site_indices() const { return site_indicies_; }

  /// Return the indices of the sites involved in the bond within a particle.
  int site(const int index) const { return site_indicies_[index]; }

  /// Return the number of sites in bond.
  int num_sites() const { return static_cast<int>(site_indicies_.size()); }

  /// Add site index.
  void add_site_index(const int index) { site_indicies_.push_back(index); }

  /// Return the name fo the bond.
  std::string name() { return name_; }

  /// Set the name of the bond.
  void set_name(const std::string name) { name_ = name; }

  void serialize(std::ostream& ostr) const;
  Bond(std::istream& istr);
  virtual ~Bond() {}

 private:
  std::string name_;
  std::vector<int> site_indicies_;
};

class Angle : public Bond {
 public:
  Angle() { set_name("angle"); }
  Angle(std::istream& istr) : Bond(istr) {}
};

class Dihedral : public Bond {
 public:
  Dihedral() { set_name("dihedral"); }
  Dihedral(std::istream& istr) : Bond(istr) {}
};

class Improper : public Bond {
 public:
  Improper() { set_name("improper"); }
  Improper(std::istream& istr) : Bond(istr) {}
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_BOND_H_
