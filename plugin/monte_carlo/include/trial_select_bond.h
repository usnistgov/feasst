
#ifndef FEASST_MONTE_CARLO_TRIAL_SELECT_BOND_H_
#define FEASST_MONTE_CARLO_TRIAL_SELECT_BOND_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select.h"

namespace feasst {

/**
  A random particle of given type is selected if previously perturbed sites
  are not available.
  Select a single bond from given anchor to mobile sites.
 */
class TrialSelectBond : public TrialSelect {
 public:
  /**
    args:
    - mobile_site: index of the mobile site.
    - anchor_site: index of the anchor site.
    - ignore_bond: if true, do not try to find a bond between mobile and anchor.
      This is used for generic reptations (default: false).
   */
  explicit TrialSelectBond(argtype args = argtype());
  explicit TrialSelectBond(argtype * arg);

  /// Return the anchor site.
  int anchor_site() const { return anchor_site_; }

  /// Return the mobile site.
  int mobile_site() const { return mobile_site_; }

  /// bond_type is added as a property.
  /// Mobile and anchor are sized.
  void precompute(System * system) override;

  bool select(const Select& perturbed,
              System * system,
              Random * random) override;

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectBond(std::istream& istr);
  virtual ~TrialSelectBond() {}

 protected:
  void serialize_trial_select_bond_(std::ostream& ostr) const;

 private:
  int mobile_site_;
  int anchor_site_;
  bool ignore_bond_;
};

inline std::shared_ptr<TrialSelectBond> MakeTrialSelectBond(
    argtype args = argtype()) {
  return std::make_shared<TrialSelectBond>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_SELECT_BOND_H_
