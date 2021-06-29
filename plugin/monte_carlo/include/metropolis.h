
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
  Metropolis();

  /// Same as above, but with an added constraint.
  Metropolis(std::shared_ptr<Constraint> constraint);

  bool is_accepted(const Acceptance& acceptance,
    const System& system,
    Random * random) override;

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<Metropolis>(istr); }
  std::shared_ptr<Criteria> create(argtype * args) const override {
    return std::make_shared<Metropolis>(); }
  void serialize(std::ostream& ostr) const override;
  Metropolis(std::istream& istr);
  ~Metropolis() {}
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
