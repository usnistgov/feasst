
#ifndef FEASST_BETA_EXPANDED_SELECT_NOTHING_H_
#define FEASST_BETA_EXPANDED_SELECT_NOTHING_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select.h"

namespace feasst {

/**
  Select nothing.
  For use with Perturb that doesn't utilize a selection.
 */
class SelectNothing : public TrialSelect {
 public:
  explicit SelectNothing(argtype args = argtype());
  explicit SelectNothing(argtype * args);

  bool select(const Select& perturbed,
              System * system,
              Random * random) override { return true; }

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit SelectNothing(std::istream& istr);
  virtual ~SelectNothing() {}
};

inline std::shared_ptr<SelectNothing> MakeSelectNothing(
    const argtype &args = argtype()) {
  return std::make_shared<SelectNothing>(args);
}

}  // namespace feasst

#endif  // FEASST_BETA_EXPANDED_SELECT_NOTHING_H_
