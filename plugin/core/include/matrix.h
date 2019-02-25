
#ifndef FEASST_CORE_MATRIX_H_
#define FEASST_CORE_MATRIX_H_

#include <vector>
#include "core/include/position.h"

namespace feasst {

/**
  A Matrix is represented by rows and columns.
  The first index is the row, and the second is the column.
 */
class Matrix {
 public:
  /// Set the size by number of rows and columns
  void set_size(const int num_rows, const int num_columns);

  /// Set the value of an element given by the row and column index.
  void set_value(const int row, const int column, const double value);

  /// Return the value of the element.
  double value(const int row, const int column) const {
    return matrix_[row][column]; }

  /// Return the entire matrix.
  std::vector<std::vector<double> > matrix() const { return matrix_; }

  /// Check that each row and column are the same size.
  virtual void check() const;

  /// Switch the rows and columns.
  void transpose();

  /// Multiply all elements by a constant.
  void multiply(const double constant);

  /// Return the vector which is a product of the multiplication of a matrix
  /// with the given vector.
  Position multiply(const Position& vec) const;

  /// Return true if all elements are equal to the given matrix.
  bool is_equal(const Matrix& matrix) const;

  /// Return the matrix as a human readable string.
  std::string str() const;

  virtual ~Matrix() {}

 private:
  std::vector<std::vector<double> > matrix_;
};

/**
  This class assumes the matrix is square with three columns and three rows.
  This assumption allows for more easy or efficient operations.
 */
class MatrixThreeByThree : public Matrix {
 public:
  /// Return the determinant.
  double determinant() const;

  /// Invert the matrix.
  void invert();

  /// Check that the matrix is square and of the assumed size.
  void check() const override;

  virtual ~MatrixThreeByThree() { check(); }
};

/**
  Rotation matrices represent a rotation of a rigid body.
  They must be square with a unit determinant.
 */
class RotationMatrix : public MatrixThreeByThree {
 // Note to HWH: extend this class for 2D as well.
 public:
  /// Compute the rotation matrix given an angle of rotation about a given
  /// axis.
  void axis_angle(const Position& axis,
    /// The angle is in units of degrees.
    const double degree_angle);

  /// Check that the derminant of a rotation matrix is unity.
  void check() const override;

  virtual ~RotationMatrix() { check(); }
};

}  // namespace feasst

#endif  // FEASST_CORE_MATRIX_H_
