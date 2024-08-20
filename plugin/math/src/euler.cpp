#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/matrix.h"
#include "math/include/euler.h"

namespace feasst {

void Euler::compute_rotation_matrix(RotationMatrix * matrix) const {
  if (matrix->num_rows() == 0) {
    matrix->set_size(3, 3);
  }
  const double s1 = std::sin(phi_);
  const double c1 = std::cos(phi_);
  const double s2 = std::sin(theta_);
  const double c2 = std::cos(theta_);
  const double s3 = std::sin(psi_);
  const double c3 = std::cos(psi_);
  matrix->set_value(0, 0, c1*c3 - c2*s1*s3);   // r11
  matrix->set_value(0, 1, -c1*s3 - c2*c3*s1);  // r12
  matrix->set_value(0, 2, s1*s2);              // r13
  matrix->set_value(1, 0, c3*s1 + c1*c2*s3);   // r21
  matrix->set_value(1, 1, -s1*s3 + c1*c2*c3);  // r22
  matrix->set_value(1, 2, -c1*s2);             // r23
  matrix->set_value(2, 0, s2*s3);              // r31
  matrix->set_value(2, 1, c3*s2);              // r32
  matrix->set_value(2, 2, c2);                 // r33
}

void Euler::set(const Matrix& matrix) {
  ASSERT(matrix.num_rows() == 3, "3D only");
  ASSERT(matrix.num_columns() == 3, "3D only");
  const std::vector<std::vector<double> >& mat = matrix.matrix();
  if (std::abs(mat[2][2] - 1.) < 1e-12) {
    theta_ = std::acos(1.);
  } else {
    theta_ = std::acos(mat[2][2]);
  }
  const double stheta = std::sin(theta_);
  if (std::abs(stheta) < 1e-8) {
    phi_ = 0.;
    psi_ = std::atan2(mat[1][0], mat[0][0]);
  } else {
    phi_ = std::atan2(mat[0][2]/stheta, -mat[1][2]/stheta);
    psi_ = std::atan2(mat[2][0]/stheta,  mat[2][1]/stheta);
  }
}

bool Euler::is_equal(const Euler& euler, const double tolerance) const {
  if (std::abs(phi_ - euler.phi()) > tolerance) {
    return false;
  } else if (std::abs(theta_ - euler.theta()) > tolerance) {
    return false;
  } else if (std::abs(psi_ - euler.psi()) > tolerance) {
    return false;
  }
  return true;
}

std::string Euler::str() const {
  std::stringstream ss;
  ss << phi_ << "," << theta_ << "," << psi_;
  return ss.str();
}

void Euler::serialize(std::ostream& ostr) const {
  feasst_serialize_version(3509, ostr);
  feasst_serialize(phi_, ostr);
  feasst_serialize(theta_, ostr);
  feasst_serialize(psi_, ostr);
}

Euler::Euler(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3509, "unrecognized version: " << version);
  feasst_deserialize(&phi_, istr);
  feasst_deserialize(&theta_, istr);
  feasst_deserialize(&psi_, istr);
}

}  // namespace feasst
