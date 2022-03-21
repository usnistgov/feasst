
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
  AlwaysReject();

  /// Same as above, but with an added constraint.
  explicit AlwaysReject(std::shared_ptr<Constraint> constraint);

  bool is_accepted(
    const System& system,
    Acceptance * acceptance,
    Random * random) override { return false; }

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<AlwaysReject>(istr); }

  void serialize(std::ostream& ostr) const override;
  AlwaysReject(std::istream& istr);
  ~AlwaysReject() {}
};

inline std::shared_ptr<AlwaysReject> MakeAlwaysReject() {
  return std::make_shared<AlwaysReject>();
}

inline std::shared_ptr<AlwaysReject> MakeAlwaysReject(
    std::shared_ptr<Constraint> constraint) {
  return std::make_shared<AlwaysReject>(constraint);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_ALWAYS_REJECT_H_
