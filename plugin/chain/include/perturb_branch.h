
#ifndef FEASST_CHAIN_PERTURB_BRANCH_H_
#define FEASST_CHAIN_PERTURB_BRANCH_H_

#include <memory>
#include <iostream>
#include "monte_carlo/include/perturb_distance_angle.h"
#include "system/include/rigid_angle.h"
#include "chain/include/perturb_branch.h"

namespace feasst {

/**
  Branching is a complex operation because of the correlation between
  multiple angles.
  See BondThreeBody for more information.
 */
class PerturbBranch : public PerturbMove {
 public:
  //@{
  /** @name Arguments
   */
  explicit PerturbBranch(argtype args = argtype());
  explicit PerturbBranch(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Same as PerturbDistanceAngle, but for two sets of angles/bonds.
  /// Also, obtain the angle parameters between the two mobile branch sites.
  void precompute(TrialSelect * select, System * system) override;

  /// Place second mobile site in branch as described in SelectBranch.
  void place_in_branch(const double distance,  // bond length
    const double angle,  // angle as in PerturbDistanceAngle::place_in_circle
    const double branch_angle,  // angle formed between both mobile sites.
    /// Angle formed by the first mobile site and both anchors.
    /// Used to check for a planar branch.
    const double extra_angle,
    System * system,
    TrialSelect * select,
    Random * random);

  void move(const bool is_position_held,
    System * system,
    TrialSelect * select,
    Random * random,
    Acceptance * acceptance) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbBranch(std::istream& istr);
  virtual ~PerturbBranch() {}

  //@}
 private:
  PerturbDistanceAngle a2a1m1_, a2a1m2_, m1a1m2_;

  /// solve equations for branch, returning particle 3 given 1 and 2
  void solve_branch_(const double x1, const double y1, const double z1,
                     const double x2, const double y2, const double z2,
                     double *x3, double *y3, double *z3, const double c143,
                     const double c243, const bool planar,
                     Random * random) const;

  // temporary
  RigidAngle angle_;
};

inline std::shared_ptr<PerturbBranch> MakePerturbBranch(
    argtype args = argtype()) {
  return std::make_shared<PerturbBranch>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_BRANCH_H_
