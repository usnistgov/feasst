
#ifndef FEASST_CONFIGURATION_CELLS_H_
#define FEASST_CONFIGURATION_CELLS_H_

#include <string>
#include <vector>
#include "configuration/include/select.h"

namespace feasst {

/**
  Divide a cuboid domain into cells.
 */
class Cells {
 public:
  Cells() { set_group(); }

  // HWH: better optimize method of building list of neighboring cells.
  // HWH: currently too slow to NPT, frequent rebuilds or large systems.
  // HWH: consider elongated boxes for minimal requirement.
  /// Create the number, length and neighbors.
  /// By default, abort if there aren't more than \f$3^D\f$ cells,
  /// where D is the dimension.
  void create(const double min_length, const std::vector<double> side_lengths);

  /// Return the number.
  int num_total() const;

  /// Return the number in a dimension.
  int num(const int dimension) const { return num_[dimension]; }

  /// Return the number.
  std::vector<int> num() const { return num_; }

  /// Set the group.
  void set_group(const int index = 0) { group_ = index; }

  /// Return the group.
  int group() const { return group_; }

  /// Clear all private member data.
  void clear();

  /// Return the neighbors.
  /// The first index is the cell.
  /// The second is a list of neighboring cells (including self).
  const std::vector<std::vector<int> >& neighbor() const { return neighbor_; }

  /// Return the particles and sites within the cells.
  /// The first index is the cell index.
  const std::vector<Select>& particles() const { return particles_; }

  /// Return the number of particles within the cells.
  int num_sites() const;

  /// Return the unique number cell in which the scaled coordinate resides.
  /// Scaled coordinates are positions divided by the respective domain size.
  int id(const std::vector<double>& scaled_coord) const;

  /// Return the type.
  int type() const { return type_; }

  /// Set the type.
  void set_type(const int type) { type_ = type; }

  ///
  void add(const Select& select, const int cell) {
    // INFO("before " << particles_[cell].str());
    particles_[cell].add(select);
    // INFO("after " << particles_[cell].str());
  }
  void remove(const Select& select, const int cell) {
    particles_[cell].remove(select);
  }
  void update(const Select& select, const int cell_new, const int cell_old) {
    particles_[cell_old].remove(select);
    particles_[cell_new].add(select);
  }

  std::string str() const;

  void serialize(std::ostream& ostr) const;
  explicit Cells(std::istream& istr);
  virtual ~Cells() {}

 private:
  int group_;
  int type_ = -1;

  // per dimension vectors
  std::vector<int> num_;

  // per cell vectors
  std::vector<std::vector<int> > neighbor_;
  std::vector<Select> particles_;  // particles for each cell

  /// Return the unique cell id number for a given cell vector.
  int id_(std::vector<int> position);

  /// Build neighbors. HWH optimize this
  void build_neighbors_2D_();
  void build_neighbors_3D_();

  /// Build list of particles in cells.
  void build_particles_();
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_CELLS_H_
