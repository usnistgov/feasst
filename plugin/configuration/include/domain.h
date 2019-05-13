
#ifndef FEASST_CONFIGURATION_DOMAIN_H_
#define FEASST_CONFIGURATION_DOMAIN_H_

#include <math.h>
#include <vector>
#include "math/include/position.h"
#include "math/include/random.h"
#include "configuration/include/cells.h"

namespace feasst {

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

  For the first (i.e., "x"), \f$ \vec{l_x} = {l_x, 0, 0}   \f$
  For the second (i.e., "y"), \f$ vec{l_y} = {xy, l_y, 0}  \f$
  For the third (i.e., "z"), \f$ \vec{l_z} = {xz, yz, l_z} \f$

  Thus, the angle, \f$\alpha\f$, between the "x" and "y" vectors is given by
  \f$ |l_x| |l_y| \cos\alpha = \vec{l_x} \cdot \vec{l_y}\f$.

  On the inner workings of Monte Carlo codes
  https://doi.org/10.1080/08927022.2013.819102
*/
class Domain {
 public:
  Domain();

  /// Get the side lengths.
  Position side_length() const { return side_length_; }

  /// Get the side length.
  double side_length(const int dimension) const {
    return side_length_.coord(dimension);
  }

  /// Set the side lengths.
  void set_side_length(const Position& side_length);

  /// Set the side length.
  void set_side_length(const int dimension, const double length) {
    side_length_.set_coord(dimension, length);
  }

  /// Set the cubic box length.
  /// Return self for chain setting.
  Domain& set_cubic(const double box_length);

  /// Set the xy tilt factor. By default it is zero.
  Domain& set_xy(const double xy = 0.);
  double xy() const { return xy_; }

  /// Set the xz tilt factor.
  Domain& set_xz(const double xz = 0.);
  double xz() const { return xz_; }

  /// Set the yz tilt factor.
  Domain& set_yz(const double yz = 0.);
  double yz() const { return yz_; }

  /// Disable periodicity in a particular dimension.
  void disable(const int dimension) { periodic_[dimension] = false; }
  bool periodic(const int dimension) const { return periodic_[dimension]; }

  /// Return the dimensionality.
  int dimension() const { return side_length_.size(); }

  /// Return the volume.
  double volume() const;

  /// Return the shift necessary to wrap the position.
  virtual Position shift(const Position& position) const;

  /// Wrap the input position.
  void wrap(Position * position) const;

  /// Return random position within boundaries.
  /// Require random number generator to keep method constant and prevent
  /// seeding to the same value over and over again.
  Position random_position(Random * random) const;

  /// Return the minimal side length
  double min_side_length() const;

  /// HWH implement check

  /// Initialize the cells according to the minimum side length.
  void init_cells(const double min_length,
    /// A group index corresponds to a group defined in configuration.
    /// If this index is set to 0 (default) use all particles and sites.
    const int group_index = 0);

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


  /// Return the unique number cell for the position.
  int cell_id(const Position& position, const Cells& cells) const;

  bool is_tilted() const { return is_tilted_; }

  /// Return true if cell list is enabled.
  bool is_cell_enabled() const {
    if (cells().size() > 0) {
      return true;
    }
    return false;
  }

  // HWH priority 1 remove r2
  // Optimized domain wrap for use in inner pair loops.
  // Not for typical users.
  void wrap_opt(const Position& pos1,
      const Position& pos2,
      Position * rel,
      double * r2) const {
    const int dimen = pos1.dimension();
    *r2 = 0;
    const std::vector<double>& side = side_length_.coord();
    std::vector<double>* dxv = (*rel).get_coord();
    if (is_tilted_) {
      wrap_triclinic_opt(pos1, pos2, rel, r2);
      return;
    } else {
      for (int dim = 0; dim < dimen; ++dim) {
        (*dxv)[dim] = pos1.coord()[dim] - pos2.coord()[dim];
        if (dim < static_cast<int>(side.size())) {
          const double side_length = side[dim];
          if (periodic_[dim]) {
            (*dxv)[dim] -= side_length*rint((*dxv)[dim]/side_length);
          }
        }
        *r2 += (*dxv)[dim]*(*dxv)[dim];
      }
    }
  }
  void wrap_triclinic_opt(const Position& pos1,
      const Position& pos2,
      Position * rel,
      double * r2) const {
    *r2 = 0;
    const std::vector<double>& side = side_length_.coord();
    std::vector<double>* dxv = (*rel).get_coord();
    if (pos1.dimension() >= 3) {
      if (2 < static_cast<int>(side.size())) {
        if (periodic_[2]) {
          const double side_length = side[2];
          const int num_wrap = rint((*dxv)[2]/side_length);
          (*dxv)[2] -= num_wrap*side_length;
          (*dxv)[1] -= num_wrap*yz_;
          (*dxv)[0] -= num_wrap*xz_;
        }
      }
      *r2 += (*dxv)[2]*(*dxv)[2];
    }
    if (1 < static_cast<int>(side.size())) {
      if (periodic_[1]) {
        const double side_length = side[1];
        const int num_wrap = rint((*dxv)[1]/side_length);
        (*dxv)[1] -= num_wrap*side_length;
        (*dxv)[0] -= num_wrap*xy_;
      }
    }
    if (0 < static_cast<int>(side.size())) {
      if (periodic_[0]) {
        const double side_length = side[0];
        (*dxv)[0] -= side_length*rint((*dxv)[0]/side_length);
      }
    }
    *r2 += (*dxv)[0]*(*dxv)[0] + (*dxv)[1]*(*dxv)[1];
  }

  void serialize(std::ostream& ostr) const;
  Domain(std::istream& istr);
  virtual ~Domain() {}

 protected:
  Position side_length_;
  double xy_, xz_, yz_;
  bool is_tilted_ = false;
  std::vector<bool> periodic_;

  /// Cell lists
  std::vector<Cells> cells_;
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_DOMAIN_H_
