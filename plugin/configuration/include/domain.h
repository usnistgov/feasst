
#ifndef FEASST_CONFIGURATION_DOMAIN_H_
#define FEASST_CONFIGURATION_DOMAIN_H_

#include <memory>
#include <vector>
#include "math/include/position.h"
#include "configuration/include/cells.h"

namespace feasst {

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
  The cuboid domain may be faster and cells are specific to cuboid domain.

  A triclinic-shaped domain is similar to cuboid-shaped except that the angles
  between the boundaries may different from ninety degrees.
  In two-dimensions, this reduces to a parallelogram.

  The triclinic periodic cell is defined by a vector for each dimension.
  This implementation is only valid for the following two- or three-dimensions.

  For the first  (i.e., "x"):

  \f$ \vec{l_x} = {l_x, 0, 0}   \f$

  For the second (i.e., "y"):

  \f$ \vec{l_y} = {xy, l_y, 0}  \f$

  For the third  (i.e., "z"):

  \f$ \vec{l_z} = {xz, yz, l_z} \f$

  Thus, the angle, \f$\alpha\f$, between the "x" and "y" vectors is given by

  \f$ |l_x| |l_y| \cos\alpha = \vec{l_x} \cdot \vec{l_y}\f$.

  On the inner workings of Monte Carlo codes
  https://doi.org/10.1080/08927022.2013.819102
*/
class Domain {
 public:
  /**
    args:
    - side_length[i]: set the side length of the i-th dimension.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
    - cubic_box_length: side length of cubic perioidic boundary conditions.
    - xy: set the tilt along the x-y direction (default: 0).
    - xz: set the tilt along the x-z direction (default: 0).
    - yz: set the tilt along the y-z direction (default: 0).
    - init_cells[i]: build cell list with given minimum length between cells.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
      If only one cell, the "[i]" is optional.
    - cell_group[i]: only compute cells for those in given group index
      for cell of corresponding "i" (default: 0).
   */
  explicit Domain(const argtype& args = argtype());

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

  /// Disable periodicity in a particular dimension.
  void disable(const int dimension) { periodic_[dimension] = false; }
  bool periodic(const int dimension) const { return periodic_[dimension]; }

  /// Return the dimensionality.
  int dimension() const { return side_lengths_.size(); }

  /// Return the volume.
  double volume() const;

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

  /// HWH implement check

  /// Initialize the cells according to the minimum side length.
  void init_cells(const double min_length,
    /// A group index corresponds to a group defined in configuration.
    /// If this index is set to 0 (default) use all particles and sites.
    const int group_index = 0);

  /// Return the number of cell lists.
  int num_cells() const { return static_cast<int>(cells_.size()); }

  /// Return the cells.
  const std::vector<Cells>& cells() const { return cells_; }

  /// Return the cells by index.
  const Cells& cells(const int index) const;

  /// Add selection to cells.
  void add_to_cell_list(const int cell_index,
                        const Select& select,
                        const int particle_cell) {
    cells_[cell_index].add(select, particle_cell);
  }

  /// Update selection in cells.
  void update_cell_list(const int cell_index,
                        const Select& select,
                        const int cell_new,
                        const int cell_old) {
    cells_[cell_index].update(select, cell_new, cell_old);
  }

  /// Remove selection from cells.
  void remove_from_cell_list(const int cell_index,
                             const Select& select,
                             const int cell) {
    cells_[cell_index].remove(select, cell);
  }


  /// Return the unique cell number for the position.
  int cell_id(const Position& position, const Cells& cells);

  bool is_tilted() const { return is_tilted_; }

  /// Return true if cell list is enabled.
  bool is_cell_enabled() const {
    if (cells().size() > 0) {
      return true;
    }
    return false;
  }

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
  std::string status_header() const;

  /// Return the brief status for periodic output.
  std::string status() const;

  void serialize(std::ostream& ostr) const;
  explicit Domain(std::istream& istr);
  virtual ~Domain() {}

 protected:
  Position side_lengths_;
  double xy_, xz_, yz_;
  bool is_tilted_ = false;
  std::vector<bool> periodic_;

  /// used for optimized pbc
  Position opt_origin_, opt_rel_, opt_pbc_;
  double opt_r2_;
  void resize_opt_(const int dimension);

  /// Cell lists
  std::vector<Cells> cells_;

  void set_xy_(const double yz);
  void set_xz_(const double yz);
  void set_yz_(const double yz);
};

inline std::shared_ptr<Domain> MakeDomain(const argtype &args = argtype()) {
  return std::make_shared<Domain>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_DOMAIN_H_
