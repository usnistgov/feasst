
#ifndef FEASST_MONTE_CARLO_TRIAL_SELECT_DIHEDRAL_H_
#define FEASST_MONTE_CARLO_TRIAL_SELECT_DIHEDRAL_H_

#include <map>
#include <string>
#include <memory>
#include "monte_carlo/include/trial_select_angle.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  A random particle of given type is selected if previously perturbed sites are
  not available.
  Select a single dihedral from three given anchor sites and one mobile site.
  The mobile site is directly bonded to the first anchor site, the second anchor
  site is bonded to first anchor, and the third anchor is bonded to the second.

  anchor3("l") - anchor2("k") - anchor("j") - mobile("i")

  See Position::torsion_angle_radians for more information.
 */
class TrialSelectDihedral : public TrialSelectAngle {
 public:
  //@{
  /** @name Arguments
    - anchor_site3 : name of third anchor site.
   */
  explicit TrialSelectDihedral(argtype args = argtype());
  explicit TrialSelectDihedral(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Same as TrialSelectAngle, but also add the third anchor site, and add
  /// dihedral_type as an anchor property.
  void precompute(System * system) override;

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectDihedral(std::istream& istr);
  virtual ~TrialSelectDihedral() {}

  //@}
 protected:
  void serialize_trial_select_dihedral_(std::ostream& ostr) const;

 private:
  std::string anchor_site3_name_;
};

inline std::shared_ptr<TrialSelectDihedral> MakeTrialSelectDihedral(
    argtype args = argtype()) {
  return std::make_shared<TrialSelectDihedral>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_SELECT_DIHEDRAL_H_
