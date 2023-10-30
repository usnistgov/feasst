
#ifndef FEASST_MONTE_CARLO_TRIAL_SELECT_ANGLE_H_
#define FEASST_MONTE_CARLO_TRIAL_SELECT_ANGLE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select_bond.h"

namespace feasst {

/**
  A random particle of given type is selected if previously perturbed sites are
  not available.
  Select a single angle from two given anchor sites and one mobile site.
  The mobile site is directly bonded to the first anchor site, and the second
  anchor site is bonded to the first anchor site.
  The angle is defined as: anchor2 - anchor1 - mobile, where anchor1 is the
  vertex.

  In 2D, angle i-j-k is not the same as angle k-j-i.
  Thus, TrialSelectAngle sites in the order mobile-anchor-anchor2 must be in the
  same order as in FileParticle.
 */
class TrialSelectAngle : public TrialSelectBond {
 public:
  //@{
  /** @name Arguments
    - anchor_site2 : index of second anchor site.
    - TrialSelectBond arguments.
   */
  explicit TrialSelectAngle(argtype args = argtype());
  explicit TrialSelectAngle(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return the second anchor site.
  int anchor_site2() const { return anchor_site2_; }

  /// Same as TrialSelectBond, but also add the second anchor site, and add
  /// angle_type as an anchor property.
  void precompute(System * system) override;

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectAngle(std::istream& istr);
  virtual ~TrialSelectAngle() {}

  //@}
 protected:
  void serialize_trial_select_angle_(std::ostream& ostr) const;

 private:
  int anchor_site2_;
};

inline std::shared_ptr<TrialSelectAngle> MakeTrialSelectAngle(
    argtype args = argtype()) {
  return std::make_shared<TrialSelectAngle>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_SELECT_ANGLE_H_
