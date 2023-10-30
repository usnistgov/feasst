
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
  Bond() { class_name_ = "Bond"; }

  /// Return the indices of the sites involved in the bond within a particle.
  const std::vector<int>& site_indices() const { return site_indicies_; }

  /// Return the indices of the sites involved in the bond within a particle.
  int site(const int index) const { return site_indicies_[index]; }

  /// Return the number of sites in bond.
  int num_sites() const { return static_cast<int>(site_indicies_.size()); }

  /// Add site index.
  void add_site_index(const int index) { site_indicies_.push_back(index); }

  /// Set the bond model.
  void set_model(const std::string model) { model_ = model; }

  /// Return the bond model.
  std::string model() const { return model_; }

  std::string class_name() { return class_name_; }

  void serialize(std::ostream& ostr) const;
  explicit Bond(std::istream& istr);
  virtual ~Bond() {}

 protected:
  std::string class_name_;

 private:
  std::vector<int> site_indicies_;
  std::string model_;
};

/**
  The Angle has three site indices, listed in order of the angle.
  For example, angle ABC has a vertex at B.
  The 0 index is A, 1 is B and 2 is C.
 */
class Angle : public Bond {
 public:
  Angle() { class_name_ = "Angle"; }
  explicit Angle(std::istream& istr) : Bond(istr) {}
};

/**
  A Dihedral has four site indicies, listed in order of the angle.
  For example, dihedral ABCD has a vertex along the BC vector.
  The 0 index is A, 1 is B, 2 is C and 3 is D.
 */
class Dihedral : public Bond {
 public:
  Dihedral() { class_name_ = "Dihedral"; }
  explicit Dihedral(std::istream& istr) : Bond(istr) {}
};

/**
  Impropers are not fully implemented.
 */
class Improper : public Bond {
 public:
  Improper() { class_name_ = "Improper"; }
  explicit Improper(std::istream& istr) : Bond(istr) {}
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_BOND_H_
