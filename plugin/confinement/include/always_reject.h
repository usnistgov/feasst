
#ifndef FEASST_CONFINEMENT_ALWAYS_REJECT_H_
#define FEASST_CONFINEMENT_ALWAYS_REJECT_H_

#include "monte_carlo/include/criteria.h"

namespace feasst {

class Random;

/**
  Always reject the trial.
 */
class AlwaysReject : public Criteria {
 public:
  explicit AlwaysReject(const argtype &args = argtype());

  /// Same as above, but with an added constraint.
  AlwaysReject(std::shared_ptr<Constraint> constraint,
    const argtype& args = argtype());

  bool is_accepted(const Acceptance& acceptance,
    const System& system,
    Random * random) override { return false; }

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<AlwaysReject>(istr); }

  void serialize(std::ostream& ostr) const override;
  AlwaysReject(std::istream& istr);
  ~AlwaysReject() {}
};

inline std::shared_ptr<AlwaysReject> MakeAlwaysReject(const argtype &args = argtype()) {
  return std::make_shared<AlwaysReject>(args);
}

inline std::shared_ptr<AlwaysReject> MakeAlwaysReject(
    std::shared_ptr<Constraint> constraint,
    const argtype &args = argtype()) {
  return std::make_shared<AlwaysReject>(constraint, args);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_ALWAYS_REJECT_H_
