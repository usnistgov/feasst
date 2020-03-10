
#ifndef FEASST_MONTE_CARLO_TRIAL_ADD_MULTIPLE_H_
#define FEASST_MONTE_CARLO_TRIAL_ADD_MULTIPLE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {

/// Attempt to add a pair of particles.
class TrialAddMultiple : public Trial {
 public:
  /**
    args:
    - particle_type0: the first type of particle to add.
    - particle_type1: the second type of particle to add.
   */
  explicit TrialAddMultiple(const argtype& args = argtype());
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialAddMultiple(std::istream& istr);
  virtual ~TrialAddMultiple() {}

  // parse the number of particle types.
  std::vector<int> ptypes(Arguments * args) const;

 protected:
  void serialize_trial_add_multiple_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialAddMultiple> MakeTrialAddMultiple(
    const argtype &args = argtype()) {
  return std::make_shared<TrialAddMultiple>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_ADD_MULTIPLE_H_
