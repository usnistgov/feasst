
#ifndef FEASST_CHAIN_SELECT_BRANCH_H_
#define FEASST_CHAIN_SELECT_BRANCH_H_

#include <memory>
#include "monte_carlo/include/trial_select_angle.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  A random particle of given type is selected if previously perturbed sites are
  not available.
  Select two mobile sites bonded to anchor 1, that also forms a third angle with
  anchor 2.

  Thus, the bonded topology appears as:

::

        anchor2
           |
        anchor1
        /     \
    mobile1   mobile2

  In the source code, these are often shortened as anchor1->a1, mobile2->m2, etc.
 */
class SelectBranch : public TrialSelectAngle {
 public:
  //@{
  /** @name Arguments
    - mobile_site2 : name of second mobile site.
    - TrialSelectAngle arguments.
   */
  explicit SelectBranch(argtype args = argtype());
  explicit SelectBranch(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Same as derived, but also add second mobile site.
  void precompute(System * system) override;

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit SelectBranch(std::istream& istr);
  virtual ~SelectBranch() {}

  //@}
 protected:
  void serialize_select_branch_(std::ostream& ostr) const;

 private:
  std::string mobile_site2_name_;
};

inline std::shared_ptr<SelectBranch> MakeSelectBranch(
    const argtype &args = argtype()) {
  return std::make_shared<SelectBranch>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_SELECT_BRANCH_H_
