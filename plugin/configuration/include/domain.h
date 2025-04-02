
#ifndef FEASST_CONFIGURATION_DOMAIN_H_
#define FEASST_CONFIGURATION_DOMAIN_H_

#include <map>
#include <memory>
#include <string>
#include <vector>
#include "math/include/position.h"
#include "math/include/matrix.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class Random;

/**
  A Domain represents the spatial boundaries and constraints imposed upon the
  positions of the particles.

  The origin is always located at the center of the domain.

  By default, periodicity in each dimension is enabled when the side lengths
  are set.

  A cuboid-shaped domain may have unique side lengths but the angles between
  the edges of the domain are fixed at ninety degrees.
  In two-dimensions, the cuboid reduces to a rectangle.

  A triclinic-shaped domain is similar to cuboid-shaped except that the angles
  between the boundaries may different from ninety degrees.
  In two-dimensions, this reduces to a parallelogram.

  The triclinic periodic cell is defined by a vector for each dimension.
  This implementation is only valid for two- or three-dimensions.

  For the first  (i.e., "x"):

  \f$ \vec{l}_x = \{l_x, 0, 0\}   \f$

  For the second (i.e., "y"):

  \f$ \vec{l}_y = \{xy, l_y, 0\}  \f$

  For the third  (i.e., "z"):

  \f$ \vec{l}_z = \{xz, yz, l_z\} \f$

  Thus, the angle, \f$\alpha\f$, between the "x" and "y" vectors is given by

  \f$ |l_x| |l_y| \cos\alpha = \vec{l}_x \cdot \vec{l}_y\f$.

  Scaled coordinates, \f$\vec{s}\f$, defined in the range of [0, 1] in each
  dimension, are related to the Cartesian coordinates, \f$\vec{x}\f$,
  via the matrix, \f$ \mathbf{H}=[\vec{l}_x,\vec{l}_y,\vec{l}_z] \f$, as follows:

  \f$ \vec{x} = \mathbf{H} \vec{s} \f$

  \f$ \vec{s} = \mathbf{H}^{-1} \vec{x} \f$

  \f$ \mathbf{H} = \begin{bmatrix} lx & xy & xz\\ 0 & ly & yz\\ 0 & 0 & lz \end{bmatrix} \f$

  \f$ \mathbf{H}^{-1} = \begin{bmatrix} 1/lx & -\frac{xy}{lx* ly} & \frac{xy* yz - xz* ly}{lx* ly* lz}\\ 0 & 1/ly & -\frac{yz}{ly* lz}\\ 0 & 0 & 1/lz \end{bmatrix} \f$

  For more information, see LAMMPS:
  https://docs.lammps.org/Howto_triclinic.html

  Another good resource is Section 4.1 of
  https://doi.org/10.1080/08927022.2013.819102

  Periodic boundaries conditions are applied in scalar coordinates as

  \f$ \vec{s}' = \vec{s} - \mathrm{rint}(\vec{s}) \f$

  where rint rounds a real number to the nearest integer.

  By comparison, the angles \f$ \alpha, \beta, \gamma \f$ and scalar boundary
  lengths \f$ a, b, c \f$ are related to the above variables as follows:

  \f$ lx = a \f$

  \f$ ly = b \cos(\gamma) \f$

  \f$ lz = c \sqrt{1-\cos(\beta)^2-\eta} \f$

  \f$ \eta = \frac{\cos(\alpha) - \cos(\gamma)\cos(\beta)}{\sin(\gamma)} \f$

  \f$ xy = ly \cos(\gamma) \f$

  \f$ xz = c \cos(\beta) \f$

  \f$ yz = c \eta \f$

  Or in the other direction

  \f$ a = lx\f$

  \f$ b^2 = ly^2 + xy^2 \f$

  \f$ c^2 = lz^2 + xz^2 + yz^2 \f$

  \f$ \cos(\alpha) = \frac{xy*xz + ly*yz}{b*c}\f$

  \f$ \cos(\beta) = xz/c\f$

  \f$ \cos(\gamma) = xy/b\f$

*/
class Domain {
 public:
  //@{
  /** @name Arguments
    - side_length[i]: set the side length of the i-th dimension.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
    - cubic_side_length: side length of cubic perioidic boundary conditions.
    - periodic[i]: set if the i-th dimension is periodic (default: true).
    - xy: set the tilt along the x-y direction (default: 0).
    - xz: set the tilt along the x-z direction (default: 0).
    - yz: set the tilt along the y-z direction (default: 0).
   */
  explicit Domain(argtype args = argtype());
  explicit Domain(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Get the side lengths.
  const Position& side_lengths() const { return side_lengths_; }

  /// Get the side length.
  double side_length(const int dimension) const {
    return side_lengths_.coord(dimension);
  }

  /// Set the side lengths.
  void set_side_lengths(const Position& side_lengths);

  /// Set the side length.
  void set_side_length(const int dimension, const double length);

  /// Add a side length (and dimensionality) to the domain.
  void add_side_length(const double length);

  /// Set the cubic box length.
  /// Return self for chain setting.
  Domain& set_cubic(const double box_length);

  /// Return true if all side lengths are equal.
  bool is_cubic() const;

  /// Return the xy tilt factor. By default it is zero.
  double xy() const { return xy_; }

  /// Return the xz tilt factor.
  double xz() const { return xz_; }

  /// Return the yz tilt factor.
  double yz() const { return yz_; }

  /// Disable periodicity in a given dimension.
  void disable(const int dimension) { periodic_[dimension] = false; }

  /// Return true if given dimension is periodic.
  bool periodic(const int dimension) const { return periodic_[dimension]; }

  /// Return the dimensionality.
  int dimension() const { return side_lengths_.size(); }

  /// Return the volume.
  double volume() const;

  /// Return the matrix to multiply scaled coordinates to obtain Cartesian
  const Matrix& h() const { return h_; }

  /// Return the matrix to multiply Cartesian coordinates to obtain scaled
  const Matrix& h_inv() const { return h_inv_; }

  /// Return the scaled coordinates, s, given cartesian cordinates, x.
  void cartesian2scaled(const Position& cartesian, Position * scaled) const;

  /// Same as above, but wrap the scaled coordinates.
  void cartesian2scaled_wrap(const Position& cartesian, Position * scaled) const;

  /// Return the shift necessary to wrap the position.
  /// This is an unoptimized version.
  Position shift(const Position& position) const;

  /// Same as above, but optimized to use pre-initialized private data.
  const Position& shift_opt(const Position& position);

  /// Wrap the input position.
  /// This is an unoptimized version.
  void wrap(Position * position) const;

  /// Return random position within boundaries.
  Position random_position(Random * random) const;

  /// Optimized verison of the above, which allows update of existing position.
  void random_position(Position * position, Random * random) const;

  /// Return the minimum side length
  double min_side_length() const;

  /// Return the maximum side length
  double max_side_length() const;

  /// Return the largest possible diameter of a sphere inscribed inside domain.
  double inscribed_sphere_diameter() const;

  // HWH implement check

  bool is_tilted() const { return is_tilted_; }

  // Optimized domain wrap for use in inner pair loops.
  // Not for typical users.
  void wrap_opt(const Position& pos1,
      const Position& pos2,
      Position * rel,
      Position * pbc,
      double * r2) const;

  void wrap_triclinic_opt(const Position& pos1,
      const Position& pos2,
      Position * rel,
      Position * pbc,
      double * r2) const;

  /// Return the shift for number of wraps, num_wrap, in a given dimension, dim.
  void unwrap(const int dim, const int num_wrap, Position * shift) const;

  /// Return the header of the status for periodic output.
  std::string status_header(const std::string append = "") const;

  /// Return the brief status for periodic output.
  std::string status() const;

  void serialize(std::ostream& ostr) const;
  explicit Domain(std::istream& istr);
  virtual ~Domain() {}

  //@}
 protected:
  Position side_lengths_;
  double xy_, xz_, yz_;
  bool is_tilted_ = false;
  std::vector<bool> periodic_;
  Matrix h_, h_inv_; // Section 4.1 of https://doi.org/10.1080/08927022.2013.819102

  /// used for optimized pbc
  Position opt_origin_, opt_rel_, opt_pbc_;
  double opt_r2_;
  void resize_opt_(const int dimension);

  void set_xy_(const double yz);
  void set_xz_(const double yz);
  void set_yz_(const double yz);
  void update_h_();
};

inline std::shared_ptr<Domain> MakeDomain(argtype args = argtype()) {
  return std::make_shared<Domain>(args); }

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_DOMAIN_H_
