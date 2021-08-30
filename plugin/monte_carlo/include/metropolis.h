
#ifndef FEASST_MONTE_CARLO_METROPOLIS_H_
#define FEASST_MONTE_CARLO_METROPOLIS_H_

#include "monte_carlo/include/criteria.h"

namespace feasst {

class Random;

/**
  Metropolis acceptance criteria.

  https://doi.org/10.1063%2F1.1699114
 */
class Metropolis : public Criteria {
 public:
  /**
    args:
    - num_attempts_per_iteration: set the number of MonteCarlo trials,
      or attempts (as measured by number of calls to is_accepted) default: 1e9.
   */
  explicit Metropolis(argtype args = argtype());

  /// Same as above, but with an added constraint.
  explicit Metropolis(std::shared_ptr<Constraint> constraint);

  bool is_accepted(const Acceptance& acceptance,
    const System& system,
    Random * random) override;

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<Metropolis>(istr); }
  std::shared_ptr<Criteria> create(argtype * args) const override {
    return std::make_shared<Metropolis>(); }
  void serialize(std::ostream& ostr) const override;
  explicit Metropolis(std::istream& istr);
  ~Metropolis() {}

 private:
  int num_attempts_per_iteration_;
};

inline std::shared_ptr<Metropolis> MakeMetropolis() {
  return std::make_shared<Metropolis>();
}

inline std::shared_ptr<Metropolis> MakeMetropolis(
    std::shared_ptr<Constraint> constraint) {
  return std::make_shared<Metropolis>(constraint);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_METROPOLIS_H_
