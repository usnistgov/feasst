
#ifndef FEASST_MATH_MATRIX_H_
#define FEASST_MATH_MATRIX_H_

#include <string>
#include <vector>
#include "math/include/position.h"

namespace feasst {

/**
  A Matrix is represented by rows and columns.
  The first index is the row, and the second is the column.
 */
class Matrix {
 public:
  Matrix() {}

  /// Set the size by number of rows and columns
  void set_size(const int num_rows, const int num_columns);

  /// Alternatively, construct with 2d vector data.
  Matrix(std::vector<std::vector<double> > matrix) { matrix_ = matrix; }

  /// Return the number of rows.
  int num_rows() const { return static_cast<int>(matrix_.size()); }

  /// Return the number of columns.
  int num_columns() const { return static_cast<int>(matrix_[0].size()); }

  /// Set the value of an element given by the row and column index.
  void set_value(const int row, const int column, const double value);

  /// Return the value of the element.
  double value(const int row, const int column) const {
    return matrix_[row][column]; }

  /// Return the entire matrix.
  const std::vector<std::vector<double> >& matrix() const { return matrix_; }

  /// Check that each row and column are the same size.
  virtual void check() const;

  /// Switch the rows and columns.
  void transpose();

  /// Multiply all elements by a constant.
  void multiply(const double constant);

  /// Set all elements to zero.
  void zero() { multiply(0.); }

  /// Return the vector which is a product of the multiplication of a matrix
  /// with the given vector.
  Position multiply(const Position& vec) const;

  /// Same as above, but optimized to avoid construction of Position.
  void multiply(const Position& vec, Position * result) const;

  /// Return true if all elements are equal to the given matrix.
  bool is_equal(const Matrix& matrix) const;

  /// Return the matrix as a human readable string.
  std::string str() const;

  /// Return the determinant (only implemented for 2x2 and 3x3 matrices).
  double determinant() const;

  /// Invert the matrix (only implemented for 2x2 and 3x3 matrices).
  void invert();

  virtual ~Matrix() {}

 private:
  std::vector<std::vector<double> > matrix_;
};

/**
  Rotation matrices represent a rotation of a rigid body.
  They must be square with a unit determinant.
 */
class RotationMatrix : public Matrix {
 public:
  /// Compute the rotation matrix given an angle of rotation about a given
  /// axis.
  /// In 2D, the axis is simply used to specify the number of dimensions.
  RotationMatrix& axis_angle(const Position& axis,
    /// The angle is in units of degrees.
    const double degree_angle);

  /// Same as above, but optimized (1) axis is assumed to be unit length and
  /// (2) rotation matrix is the correct size and (3) determinant is not
  /// checked.
  void axis_angle_opt(const Position& unit_axis,
    /// The angle is in units of degrees.
    const double degree_angle);

  /// Rotate a position given a pivot point.
  void rotate(const Position& pivot, Position * rotated) const;

  /// Check square matrix, unit derminant, in addition to Matrix::check.
  void check() const override;

  virtual ~RotationMatrix() {}
};

}  // namespace feasst

#endif  // FEASST_MATH_MATRIX_H_
