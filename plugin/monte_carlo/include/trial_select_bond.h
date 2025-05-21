
#ifndef FEASST_MONTE_CARLO_TRIAL_SELECT_BOND_H_
#define FEASST_MONTE_CARLO_TRIAL_SELECT_BOND_H_

#include <map>
#include <string>
#include <memory>
#include "monte_carlo/include/trial_select.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  A random particle of given type is selected if previously perturbed sites
  are not available.
  Select a single bond from given anchor to mobile sites.
 */
class TrialSelectBond : public TrialSelect {
 public:
  //@{
  /** @name Arguments
    - mobile_site: name of the mobile site.
    - anchor_site: name of the anchor site.
    - ignore_bond: if true, do not try to find a bond between mobile and anchor.
      This is used for generic reptations (default: false).
    - TrialSelect arguments.
   */
  explicit TrialSelectBond(argtype args = argtype());
  explicit TrialSelectBond(argtype * arg);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return the anchor site.
  const std::string& anchor_site_name() const { return anchor_site_name_; }

  /// Return the mobile site.
  const std::string& mobile_site_name() const { return mobile_site_name_; }

  /// bond_type is added as a property.
  /// Mobile and anchor are sized.
  void precompute(System * system) override;

  bool select(const Select& perturbed,
              System * system,
              Random * random,
              TrialSelect * previous_select) override;

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectBond(std::istream& istr);
  virtual ~TrialSelectBond() {}

  //@}
 protected:
  void serialize_trial_select_bond_(std::ostream& ostr) const;

 private:
  std::string mobile_site_name_, anchor_site_name_;
  bool ignore_bond_;
};

inline std::shared_ptr<TrialSelectBond> MakeTrialSelectBond(
    argtype args = argtype()) {
  return std::make_shared<TrialSelectBond>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_SELECT_BOND_H_
