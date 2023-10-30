
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
  //@{
  /** @name Arguments
    - num_trials_per_iteration: define an iteration as a number of trials
      (as measured by number of calls to is_accepted) default: 1e9.
      Note that iterations are defined like cycles, but are not necessarily
      the number of particles.
    - Criteria arguments.
   */
  explicit Metropolis(argtype args = argtype());
  explicit Metropolis(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Same as above, but with an added constraint.
  explicit Metropolis(std::shared_ptr<Constraint> constraint);

  bool is_accepted(
    const System& system,
    Acceptance * acceptance,
    Random * random) override;

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<Metropolis>(istr); }
  std::shared_ptr<Criteria> create(argtype * args) const override {
    return std::make_shared<Metropolis>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Metropolis(std::istream& istr);
  ~Metropolis() {}

  //@}
 private:
  int num_trials_per_iteration_;
};

inline std::shared_ptr<Metropolis> MakeMetropolis(argtype args = argtype()) {
  return std::make_shared<Metropolis>(args);
}

inline std::shared_ptr<Metropolis> MakeMetropolis(
    std::shared_ptr<Constraint> constraint) {
  return std::make_shared<Metropolis>(constraint);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_METROPOLIS_H_
