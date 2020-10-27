
#ifndef FEASST_CHAIN_SELECT_BRANCH_H_
#define FEASST_CHAIN_SELECT_BRANCH_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select_angle.h"

namespace feasst {

/**
  A random particle of given type is selected if previously perturbed sites are
  not available.
  Select two mobile sites bonded to anchor 1, that also forms a third angle with
  anchor 2.
\rst
Thus, the bonded topology appears as:
    anchor2
       |
    anchor1
    /     \
mobile1   mobile2
\endrst
 */
class SelectBranch : public TrialSelectAngle {
 public:
  /**
    args:
    - mobile_site2 : index of second mobile site.
   */
  explicit SelectBranch(const argtype& args = argtype());

  /// Same as derived, but also add second mobile site.
  void precompute(System * system) override;

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit SelectBranch(std::istream& istr);
  virtual ~SelectBranch() {}

 protected:
  void serialize_select_branch_(std::ostream& ostr) const;

 private:
  int mobile_site2_;
};

inline std::shared_ptr<SelectBranch> MakeSelectBranch(
    const argtype &args = argtype()) {
  return std::make_shared<SelectBranch>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_SELECT_BRANCH_H_
