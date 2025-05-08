
#ifndef FEASST_MONTE_CARLO_TRIAL_SELECT_ALL_H_
#define FEASST_MONTE_CARLO_TRIAL_SELECT_ALL_H_

#include <map>
#include <string>
#include <vector>
#include <memory>
#include "monte_carlo/include/trial_select.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/// Select all particles.
class TrialSelectAll : public TrialSelect {
 public:
  explicit TrialSelectAll(argtype args = argtype());
  explicit TrialSelectAll(argtype * args);

  bool select(const Select& perturbed,
              System* system,
              Random * random,
              TrialSelect * previous_select) override;

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectAll(std::istream& istr);
  virtual ~TrialSelectAll() {}

 protected:
  void serialize_trial_select_all_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialSelectAll> MakeTrialSelectAll(
    argtype args = argtype()) {
  return std::make_shared<TrialSelectAll>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_SELECT_ALL_H_
