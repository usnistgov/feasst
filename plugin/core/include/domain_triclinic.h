
#ifndef FEASST_CORE_DOMAIN_TRICLINIC_H_
#define FEASST_CORE_DOMAIN_TRICLINIC_H_

#include "core/include/domain.h"

namespace feasst {

/**
  A triclinic-shaped domain is similar to cuboid-shaped except that the angles
  between the boundaries may different from ninety degrees.

  In this implementation, the triclinic periodic cell is defined by a vector for
  each dimension.
  This implementation is only valid for the following two- or three-dimensions.

  For the first (i.e., "x"), \f$ \vec{l_x} = {l_x, 0, 0}   \f$
  For the second (i.e., "y"), \f$ vec{l_y} = {xy, l_y, 0}  \f$
  For the third (i.e., "z"), \f$ \vec{l_z} = {xz, yz, l_z} \f$

  Thus, the angle, \f$\alpha\f$, between the "x" and "y" vectors is given by
  \f$ |l_x| |l_y| \cos\alpha = \vec{l_x} \cdot \vec{l_y}\f$.
 */
class DomainTriclinic : public Domain {
 public:
  DomainTriclinic() { set_xy().set_xz().set_yz(); }

  /// Set the xy tilt factor.
  DomainTriclinic& set_xy(const double xy = 0.) { xy_ = xy; return *this; }

  /// Set the xz tilt factor.
  DomainTriclinic& set_xz(const double xz = 0.) { xz_ = xz; return *this; }

  /// Set the yz tilt factor.
  DomainTriclinic& set_yz(const double yz = 0.) { yz_ = yz; return *this; }

  // see base class for comment
  Position shift(const Position& position) const override;

  // see base class for comment
  Position random_position(Random * random) const override;

  virtual ~DomainTriclinic() {}
 private:
  double xy_, xz_, yz_;
};

}  // namespace feasst

#endif  // FEASST_CORE_DOMAIN_TRICLINIC_H_
