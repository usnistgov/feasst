
#ifndef FEASST_CORE_DOMAIN_H_
#define FEASST_CORE_DOMAIN_H_

#include <vector>
#include "core/include/position.h"
#include "core/include/random.h"
#include "core/include/cells.h"

namespace feasst {

/**
  A Domain represents the spatial boundaries and constraints imposed upon the
  positions of the particles.

  The origin is always located at the center of the domain.

  By default, periodicity in each dimension is enabled.

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

  /// Set the xz tilt factor.
  Domain& set_xz(const double xz = 0.);

  /// Set the yz tilt factor.
  Domain& set_yz(const double yz = 0.);

  /// Disable periodicity in a particular dimension.
  void disable(const int dimension) { periodic_[dimension] = false; }

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

  /// HWH implement check_size

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
        if (periodic_[dim]) {
          if ((*dxv)[dim] >  0.5*side[dim]) {
            (*dxv)[dim] -= side[dim];
          } else if ((*dxv)[dim] < -0.5*side[dim]) {
            (*dxv)[dim] += side[dim];
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
      if (periodic_[2] && std::abs((*dxv)[2]) > 0.5*side[2]) {
        if ((*dxv)[2] < 0.) {
          (*dxv)[2] += side[2];
          (*dxv)[1] += yz_;
          (*dxv)[0] += xz_;
        } else {
          (*dxv)[2] -= side[2];
          (*dxv)[1] -= yz_;
          (*dxv)[0] -= xz_;
        }
      }
      *r2 += (*dxv)[2]*(*dxv)[2];
    }
    if (periodic_[1] && std::abs((*dxv)[1]) > 0.5*side[1]) {
      if ((*dxv)[1] < 0.) {
        (*dxv)[1] += side[1];
        (*dxv)[0] += xy_;
      } else {
        (*dxv)[1] -= side[1];
        (*dxv)[0] -= xy_;
      }
    }
    if (periodic_[0] && std::abs((*dxv)[0]) > 0.5*side[0]) {
      if ((*dxv)[0] < 0.) {
        (*dxv)[0] += side[0];
      } else {
        (*dxv)[0] -= side[0];
      }
    }
    *r2 += (*dxv)[0]*(*dxv)[0] + (*dxv)[1]*(*dxv)[1];
  }
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

#endif  // FEASST_CORE_DOMAIN_H_
