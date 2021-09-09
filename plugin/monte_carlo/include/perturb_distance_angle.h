
#ifndef FEASST_MONTE_CARLO_PERTURB_DISTANCE_ANGLE_H_
#define FEASST_MONTE_CARLO_PERTURB_DISTANCE_ANGLE_H_

#include "math/include/matrix.h"
#include "system/include/rigid_angle.h"
#include "monte_carlo/include/perturb_distance.h"

namespace feasst {

/**
  Put first site in selection, i, in a truncated conical shell that is formed
  by the intersection of a spherical shell about the first site in the anchor,
  and a conical shell with axis of symmetry along the vector connecting the
  first anchor site to the second anchor site.
  For angle potentials, the equilibrium angle and spring constant are as
  described in Random::bond_angle
  Currently implemented for harmonic bonds (exponent: 2), but could add an
  optional exponent model parameter to generalize this in the future.
 */
class PerturbDistanceAngle : public PerturbDistance {
 public:
  explicit PerturbDistanceAngle(argtype args = argtype());
  explicit PerturbDistanceAngle(argtype * args);

  /// Same as PerturbDistance, but also obtain the equilibrium angle and
  /// spring constant.
  /// If 2D, angles are positive when clockwise.
  /// Thus, when 2D and reverse (e.g., kji instead of ijk), angle = 2pi - angle.
  void precompute(TrialSelect * select, System * system) override;

  /// Return the angle.
  double angle_type() const { return angle_type_; }

  /// Return the randomly selected angle from the potential.
  double random_angle_radians(const System& system,
    const TrialSelect * select,
    Random * random);
  double random_angle_radians(const System& system,
    const TrialSelect * select,
    Random * random,
    double * bond_energy  // return the bond energy for Rosenbluth exclusion
  );

  /// Place mobile site randomly in the circle about the anchors.
  void place_in_circle(const double distance, const double angle,
    System * system,
    TrialSelect * select,
    Random * random);

  void move(System * system,
    TrialSelect * select,
    Random * random) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbDistanceAngle(std::istream& istr);
  virtual ~PerturbDistanceAngle() {}

 protected:
  void serialize_perturb_distance_angle_(std::ostream& ostr) const;

 private:
  int angle_type_ = 0.;

  // temporary
  Position rjk_;
  Position orthogonal_jk_;
  Position origin_;
  RotationMatrix rot_mat_;
  RigidAngle angle_;
};

inline std::shared_ptr<PerturbDistanceAngle> MakePerturbDistanceAngle(
    argtype args = argtype()) {
  return std::make_shared<PerturbDistanceAngle>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_DISTANCE_ANGLE_H_
