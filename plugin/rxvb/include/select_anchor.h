
#ifndef FEASST_RXNAVB_SELECT_ANCHOR_H_
#define FEASST_RXNAVB_SELECT_ANCHOR_H_

#include "monte_carlo/include/trial_select.h"

namespace feasst {

/// Select the entire particle of the anchor of the previous TrialSelect
class SelectAnchor : public TrialSelect {
 public:
  explicit SelectAnchor(argtype args = argtype());
  explicit SelectAnchor(argtype * args);
  bool select(const Select& perturbed,
    System* system,
    Random * random,
    TrialSelect * previous_select) override;
  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit SelectAnchor(std::istream& istr);
  virtual ~SelectAnchor() {}
};

}  // namespace feasst

#endif  // FEASST_RXNAVB_SELECT_ANCHOR_H_
