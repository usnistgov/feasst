
#ifndef FEASST_CORE_BOND_H_
#define FEASST_CORE_BOND_H_

#include <string>
#include <vector>
#include "core/include/typed_entity.h"
#include "core/include/properties.h"

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

  virtual ~Bond() {}

 private:
  std::vector<int> site_indicies_;
  std::string name_;
};

class Angle : public Bond {
 public:
  Angle() { set_name("angle"); }
};

class Dihedral : public Bond {
 public:
  Dihedral() { set_name("dihedral"); }
};

class Improper : public Bond {
 public:
  Improper() { set_name("improper"); }
};

}  // namespace feasst

#endif  // FEASST_CORE_BOND_H_
