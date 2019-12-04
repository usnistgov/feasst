
#ifndef FEASST_MONTE_CARLO_PERTURB_DISTANCE_ANGLE_H_
#define FEASST_MONTE_CARLO_PERTURB_DISTANCE_ANGLE_H_

#include "monte_carlo/include/perturb_distance.h"

namespace feasst {

/**
  Put first site in selection, i, in a sphere about the first site in anchor,
  j, and at an angle i,j,k (vertex: j) about the second site in anchor, j.
 */
class PerturbDistanceAngle : public PerturbDistance {
 public:
  explicit PerturbDistanceAngle(const argtype& args = argtype());

  void precompute(TrialSelect * select, System * system) override;

  void move(System * system,
      TrialSelect * select,
      Random * random) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbDistanceAngle(std::istream& istr);
  virtual ~PerturbDistanceAngle() {}

 private:
  double angle_ = 0.;

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
