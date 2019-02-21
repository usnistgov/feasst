
#ifndef FEASST_CORE_MATRIX_H_
#define FEASST_CORE_MATRIX_H_

#include "core/include/position.h"

namespace feasst {

class Matrix {
 public:
  void set_size(const int row, const int column) {
    matrix_.resize(row, std::vector<double>(column));
  }
  void set_value(const int row, const int column, const double value) {
    matrix_[row][column] = value;
  }
  double value(const int row, const int column) const { return matrix_[row][column]; }
  std::vector<std::vector<double> > matrix() const { return matrix_; }
  virtual void check() const;
  void transpose();
  void multiply(const double constant);
  Position multiply(const Position& vec) const;
  virtual ~Matrix() {}
 private:
  std::vector<std::vector<double> > matrix_;
};

class MatrixThreeByThree : public Matrix {
 public:
  double determinant() const;
  void invert();
  virtual void check() const override;
  virtual ~MatrixThreeByThree() {}
};

/// HWH extend for 2x2 rotation matrices
class RotationMatrix : public MatrixThreeByThree {
 public:
  void check() const override;
  void axis_angle(const Position& axis,
      /// angle in degrees
      const double angle);
  virtual ~RotationMatrix() {}
};

}  // namespace feasst

#endif  // FEASST_CORE_MATRIX_H_
