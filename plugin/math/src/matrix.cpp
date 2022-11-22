#include <cmath>
#include <vector>
#include "math/include/matrix.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "utils/include/io.h"

namespace feasst {

void Matrix::set_size(const int num_rows, const int num_columns) {
  matrix_.resize(num_rows, std::vector<double>(num_columns));
  // for (std::vector<double> & row : matrix_) {
  //   row.resize(num_columns);
  // }
}

void Matrix::set_value(const int row, const int column, const double value) {
  matrix_[row][column] = value;
}

void Matrix::check() const {
  ASSERT(matrix().size() != 0, "empty");
  size_t row_size = matrix()[0].size();
  for (const std::vector<double>& row : matrix_) {
    ASSERT(row.size() != 0, "empty");
    ASSERT(row.size() == row_size, "uneven");
  }
}

void Matrix::transpose(Matrix * tmp) {
  tmp->matrix_ = matrix_;
  set_size(tmp->matrix()[0].size(), tmp->num_rows());
  for (int row = 0; row < tmp->num_rows(); ++row) {
    for (int col = 0; col < tmp->num_columns(); ++col) {
      matrix_[col][row] = tmp->matrix()[row][col];
    }
  }
}

void Matrix::transpose() {
  Matrix tmp;
  transpose(&tmp);
}

void Matrix::multiply(const double constant) {
  for (int row = 0; row < static_cast<int>(matrix_.size()); ++row) {
    for (int col = 0; col < static_cast<int>(matrix_[row].size()); ++col) {
      matrix_[row][col] *= constant;
    }
  }
}

Position Matrix::multiply(const Position& vec) const {
  Position result = vec;
  multiply(vec, &result);
  return result;
}

void Matrix::multiply(const Position& vec, Position * result) const {
  for (int row = 0; row < static_cast<int>(matrix_.size()); ++row) {
    result->set_coord(row, vec.dot_product(matrix_[row]));
  }
}

bool Matrix::is_equal(const Matrix& matrix2) const {
  check();
  matrix2.check();
  if (matrix().size() != matrix2.matrix().size() ||
      matrix()[0].size() != matrix2.matrix()[0].size()) {
    return false;
  }
  for (int row = 0; row < static_cast<int>(matrix_.size()); ++row) {
    for (int col = 0; col < static_cast<int>(matrix_[row].size()); ++col) {
      if (std::abs(matrix()[row][col] -
                   matrix2.matrix()[row][col]) > NEAR_ZERO) {
        return false;
      }
    }
  }
  return true;
}

std::string Matrix::str() const { return feasst_str(matrix_); }

// thanks to https://en.wikipedia.org/wiki/Determinant
double Matrix::determinant() const {
  if (matrix_.size() == 2) {
    if (matrix_[0].size() == 2) {
      return value(0, 0)*value(1, 1) - value(0, 1)*value(1, 0);
    }
  } else if (matrix_.size() == 3) {
    if (matrix_[0].size() == 3) {
      const double a = value(0, 0);
      const double b = value(0, 1);
      const double c = value(0, 2);
      const double d = value(1, 0);
      const double e = value(1, 1);
      const double f = value(1, 2);
      const double g = value(2, 0);
      const double h = value(2, 1);
      const double i = value(2, 2);
      return a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h;
    }
  }
  FATAL("unrecognized matrix size: " << matrix_.size() << " x "
    << matrix_[0].size());
}

void Matrix::invert() {
  const double det = determinant();
  ASSERT(std::abs(det) > NEAR_ZERO, "not invertible");
  transpose();
  multiply(1./det);
}

void RotationMatrix::check() const {
  Matrix::check();
  ASSERT(matrix().size() == matrix()[0].size(), str() << "not square");
  ASSERT(std::abs(determinant() - 1.) < 10*NEAR_ZERO, str() << "not unit determinant("
    << determinant() << ")");
}

RotationMatrix& RotationMatrix::axis_angle(const Position& axis,
    const double degree_angle) {
  Position unit_axis = axis;
  if (unit_axis.size() != 2) {
    unit_axis.normalize();
  }
  set_size(unit_axis.size(), unit_axis.size());
  axis_angle_opt(unit_axis, degree_angle);
  check();
  return *this;
}

// thanks to https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
void RotationMatrix::axis_angle_opt(const Position& unit_axis,
    const double degree_angle) {
  const double radian_angle = degrees_to_radians(degree_angle);
  const double s = std::sin(radian_angle);
  const double c = std::cos(radian_angle);
  if (unit_axis.size() == 2) {
    set_value(0, 0, c);
    set_value(0, 1, -s);
    set_value(1, 0, s);
    set_value(1, 1, c);
  } else if (unit_axis.size() == 3) {
    const double C = 1-c;
    const double x = unit_axis.coord(0);
    const double y = unit_axis.coord(1);
    const double z = unit_axis.coord(2);
    set_value(0, 0, x*x*C + c);
    set_value(0, 1, x*y*C - z*s);
    set_value(0, 2, x*z*C + y*s);
    set_value(1, 0, y*x*C + z*s);
    set_value(1, 1, y*y*C + c);
    set_value(1, 2, y*z*C - x*s);
    set_value(2, 0, z*x*C - y*s);
    set_value(2, 1, z*y*C + x*s);
    set_value(2, 2, z*z*C + c);
  } else {
    FATAL("unrecognized dimension: " << unit_axis.size());
  }
}

void RotationMatrix::rotate(const Position& pivot, Position * rotated) const {
  rotated->subtract(pivot);
  *rotated = multiply(*rotated);
  rotated->add(pivot);
}

void RotationMatrix::quaternion(const Position& q) {
  ASSERT(q.size() == 4, "size: " << q.size() << 
    " but assumes 3D space, 4D qernion");
  set_value(0, 0, q.coord(0)*q.coord(0) - q.coord(1)*q.coord(1) - q.coord(2)*q.coord(2) + q.coord(3)*q.coord(3));
  set_value(1, 0, 2*(q.coord(0)*q.coord(1) + q.coord(2)*q.coord(3)));
  set_value(2, 0, 2*(q.coord(2)*q.coord(0) - q.coord(1)*q.coord(3)));
  set_value(0, 1, 2*(q.coord(0)*q.coord(1) - q.coord(2)*q.coord(3)));
  set_value(1, 1, q.coord(1)*q.coord(1) - q.coord(2)*q.coord(2) - q.coord(0)*q.coord(0) + q.coord(3)*q.coord(3));
  set_value(2, 1, 2*(q.coord(1)*q.coord(2) + q.coord(0)*q.coord(3)));
  set_value(0, 2, 2*(q.coord(2)*q.coord(0) + q.coord(1)*q.coord(3)));
  set_value(1, 2, 2*(q.coord(1)*q.coord(2) - q.coord(0)*q.coord(3)));
  set_value(2, 2, q.coord(2)*q.coord(2) - q.coord(0)*q.coord(0) - q.coord(1)*q.coord(1) + q.coord(3)*q.coord(3));
}

void Matrix::multiply(const Matrix& matrix, Matrix * result, Position * tmp1, Position * tmp2) const {
  if (result->num_rows() == 0) {
    result->set_size(num_rows(), matrix.num_columns());
  }
  if (tmp1->size() == 0) {
    tmp1->set_to_origin(num_columns());
    tmp2->set_to_origin(num_columns());
  }
  ASSERT(num_columns() == matrix.num_rows(),
    num_columns() << " != " << matrix.num_rows());
  ASSERT(result->num_rows() == num_rows(),
    result->num_rows() << " != " << num_rows());
  ASSERT(result->num_columns() == matrix.num_columns(),
    result->num_columns() << " != " << matrix.num_columns());
  Position row(num_columns());
  Position col(matrix.num_rows());
  for (int row1 = 0; row1 < num_rows(); ++row1) {
    for (int col2 = 0; col2 < matrix.num_columns(); ++col2) {
      row.set_vector(matrix_[row1]);
      for (int row2 = 0; row2 < matrix.num_rows(); ++row2) {
        col.set_coord(row2, matrix.matrix()[row2][col2]);
      }
      result->set_value(row1, col2, row.dot_product(col));
    }
  }
}

Matrix Matrix::multiply(const Matrix& matrix) const {
  Matrix result;
  Position tmp1, tmp2;
  multiply(matrix, &result, &tmp1, &tmp2);
  return result;
}

bool Matrix::is_identity(const double tolerance) const {
  for (int row = 0; row < num_rows(); ++row) {
    for (int col = 0; col < num_columns(); ++col) {
      if (row == col) {
        if (std::abs(matrix_[row][col] - 1) > tolerance) {
          return false;
        }
      } else {
        if (std::abs(matrix_[row][col]) > tolerance) {
          return false;
        }
      }
    }
  }
  return true;
}

}  // namespace feasst
