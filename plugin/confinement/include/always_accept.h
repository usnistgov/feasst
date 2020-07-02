
#ifndef FEASST_CONFINEMENT_ALWAYS_ACCEPT_H_
#define FEASST_CONFINEMENT_ALWAYS_ACCEPT_H_

#include "monte_carlo/include/criteria.h"

namespace feasst {

/**
  Always accept the trial.
 */
class AlwaysAccept : public Criteria {
 public:
  explicit AlwaysAccept(const argtype &args = argtype());

  /// Same as above, but with an added constraint.
  AlwaysAccept(std::shared_ptr<Constraint> constraint,
    const argtype& args = argtype());

  bool is_accepted(const Acceptance& acceptance,
    const System& system,
    const double uniform_random) override;

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<AlwaysAccept>(istr); }

  void serialize(std::ostream& ostr) const override;
  AlwaysAccept(std::istream& istr);
  ~AlwaysAccept() {}
};

inline std::shared_ptr<AlwaysAccept> MakeAlwaysAccept(const argtype &args = argtype()) {
  return std::make_shared<AlwaysAccept>(args);
}

inline std::shared_ptr<AlwaysAccept> MakeAlwaysAccept(
    std::shared_ptr<Constraint> constraint,
    const argtype &args = argtype()) {
  return std::make_shared<AlwaysAccept>(constraint, args);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_ALWAYS_ACCEPT_H_
