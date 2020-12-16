
#ifndef FEASST_MONTE_CARLO_PERTURB_DISTANCE_ANGLE_H_
#define FEASST_MONTE_CARLO_PERTURB_DISTANCE_ANGLE_H_

#include "math/include/matrix.h"
#include "monte_carlo/include/perturb_distance.h"

namespace feasst {

/**
 * Put first site in selection, i, in a sphere about the first site in anchor,
 * j, and at an angle i,j,k (vertex: j) about the second site in anchor, j.
 * For angle potentials, the equilibrium angle and spring constant are as
 * described in Random::bond_angle
 * Currently implemented for harmonic bonds (exponent: 2), but could add an
 * optional exponent model parameter to generalize this.
 */
class PerturbDistanceAngle : public PerturbDistance {
 public:
  explicit PerturbDistanceAngle(const argtype& args = argtype());

  /// Same as PerturbDistance, but also obtain the equilibrium angle and
  /// spring constant.
  /// If 2D, angles are positive when clockwise.
  /// Thus, when 2D and reverse (e.g., kji instead of ijk), angle = 2pi - angle.
  void precompute(TrialSelect * select, System * system) override;

  /// Return the angle.
  double angle() const { return angle_; }

  /// Return the spring constant. If rigid, return -1 (default).
  double spring_constant() const { return spring_constant_; }

  /// Return true if angle is rigid (e.g., if spring_constant is -1).
  bool is_rigid() const;

  /// Return true if angle is freely-jointed (e.g., if spring_constant is 0).
  bool is_freely_jointed() const;

  /// Return the randomly selected angle from the potential.
  /// If the spring constant is -1 (rigid), simplly return the angle.
  /// IF the spring constant is 0 (freely-jointed), use angle() for minimum.
  double random_angle(Random * random,
    const double beta,  /// inverse temperature
    const int dimension) const;

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
  double angle_ = 0.;
  double spring_constant_ = -1;

  // temporary
  Position rjk_;
  Position orthogonal_jk_;
  Position origin_;
  RotationMatrix rot_mat_;
};

inline std::shared_ptr<PerturbDistanceAngle> MakePerturbDistanceAngle(
    const argtype& args = argtype()) {
  return std::make_shared<PerturbDistanceAngle>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_DISTANCE_ANGLE_H_
