
#ifndef FEASST_CHARGE_TRIAL_ADD_MULTIPLE_H_
#define FEASST_CHARGE_TRIAL_ADD_MULTIPLE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {


// parse the number of particle types.
std::vector<int> ptypes(argtype * args);

/**
  Attempt to add multiple particles.
  Typically requires the use of a reference index.

  args:
  - particle_type[i]: the i-th type of particle to add.
    The "[i]" is to be substituted for an integer 0, 1, 2, ...
 */
class TrialAddMultiple : public Trial {
 public:
  TrialAddMultiple(argtype args = argtype());
  TrialAddMultiple(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialAddMultiple>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialAddMultiple>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialAddMultiple(std::istream& istr);
  virtual ~TrialAddMultiple() {}
};

inline std::shared_ptr<TrialAddMultiple> MakeTrialAddMultiple(argtype args = argtype()) {
  return std::make_shared<TrialAddMultiple>(args); }

}  // namespace feasst

#endif  // FEASST_CHARGE_TRIAL_ADD_MULTIPLE_H_
