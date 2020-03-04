
#ifndef FEASST_MONTE_CARLO_TRIAL_ADD_PAIR_H_
#define FEASST_MONTE_CARLO_TRIAL_ADD_PAIR_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {

/// Attempt to add a pair of particles.
class TrialAddPair : public Trial {
 public:
  /**
    args:
    - particle_type0: the first type of particle to add.
    - particle_type1: the second type of particle to add.
   */
  explicit TrialAddPair(const argtype& args = argtype());
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialAddPair(std::istream& istr);
  virtual ~TrialAddPair() {}

 protected:
  void serialize_trial_add_pair_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialAddPair> MakeTrialAddPair(
    const argtype &args = argtype()) {
  return std::make_shared<TrialAddPair>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_ADD_PAIR_H_
