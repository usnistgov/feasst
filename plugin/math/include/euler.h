
#ifndef FEASST_MATH_EULER_H_
#define FEASST_MATH_EULER_H_

#include <string>
#include <vector>

namespace feasst {

class Matrix;
class RotationMatrix;

/**
  There are many ambiguities in Euler angle and rotation matrix definitions.
  See https://en.wikipedia.org/wiki/Euler_angles
  and https://en.wikipedia.org/wiki/Rotation_matrix#Ambiguities

  FEASST uses "so-called" proper, active intrinsic z-x-z or "x-convention"
  which is given by the following three intrinsic rotations:

  1. rotate by phi [-pi, pi] about the z-axis.
  2. rotate by theta [0, pi] about the new x-axis.
  3. rotate by psi [-pi, pi] about the new z-axis.

  The rotation matrixes for each of these are as described in
  https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations

  For a z-axis rotation is given by
  [c, -s, 0
   s, c, 0
   0, 0, 1]
  where c is cosine of the angle of rotation, and s is the sine.

  and an x-axis rotation is given by
  [1, 0, 0
   0, c, -s
   0, s, c]

  The following rotation matrix is obtained for the rotation described above
  [c1c3 - c2s1s3, -c1s3 - c2c3s1, s1s2
   c3s1 + c1c2s3, c1c2c3 - s1s3, -c1s2
   s2s3, c3s2, c2]
  which should be equivalent to the Z1X2Z3 Proper Euler angle described in
  https://en.wikipedia.org/wiki/Euler_angles

  Note that these resulting rotation matrix is the inverse of the one described
  in the following website: https://mathworld.wolfram.com/EulerAngles.html.
  This is likely because the rotation matrix is a passive rotation, which
  means that it is from the perspective of changing the coordinate system
  rather than the coordinates of a particle in a fixed frame.
  https://en.wikipedia.org/wiki/Active_and_passive_transformation
 */
class Euler {
 public:
  Euler() {}

  /// Construct from values.
  Euler(const double phi, const double theta, const double psi) {
    set(phi, theta, psi); }

  /// Set the Euler angles.
  void set(const double phi, const double theta, const double psi) {
    phi_ = phi; theta_ = theta; psi_ = psi; }

  /// Set the Euler angles from RotationMatrix.
  void set(const Matrix& matrix);

  /// Return the first angle [-pi, pi] about the z-axis.
  double phi() const { return phi_; }

  /// Return the second angle [0, pi] about the new x-axis.
  double theta() const { return theta_; }

  /// Return the third angle [-pi, pi] about the new z-axis.
  double psi() const { return psi_; }

  /// Compute a RotationMatrix from Euler angles.
  void compute_rotation_matrix(RotationMatrix * matrix) const;

  /// Return true if equal.
  bool is_equal(const Euler& euler, const double tolerance = 1e-8) const;

  /// Return human readable format.
  std::string str() const;

  void serialize(std::ostream& ostr) const;
  explicit Euler(std::istream& istr);

 private:
  double phi_ = 0.;
  double theta_ = 0.;
  double psi_ = 0.;
};

}  // namespace feasst

#endif  // FEASST_MATH_EULER_H_
