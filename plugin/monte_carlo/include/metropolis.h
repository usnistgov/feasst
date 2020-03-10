
#ifndef FEASST_MONTE_CARLO_METROPOLIS_H_
#define FEASST_MONTE_CARLO_METROPOLIS_H_

#include "monte_carlo/include/criteria.h"

namespace feasst {

/**
  Metropolis acceptance criteria.
 */
class Metropolis : public Criteria {
 public:
  Metropolis(const argtype &args = argtype()) : Criteria(args) {
    class_name_ = "Metropolis";
  }

  bool is_accepted(const Acceptance& acceptance,
    const System * system,
    const double uniform_random) override;

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<Metropolis>(istr); }

  void serialize(std::ostream& ostr) const override;
  Metropolis(std::istream& istr);
  ~Metropolis() {}
};

inline std::shared_ptr<Metropolis> MakeMetropolis(const argtype &args = argtype()) {
  return std::make_shared<Metropolis>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_METROPOLIS_H_
