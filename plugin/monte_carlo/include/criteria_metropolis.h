
#ifndef FEASST_MONTE_CARLO_CRITERIA_METROPOLIS_H_
#define FEASST_MONTE_CARLO_CRITERIA_METROPOLIS_H_

#include "monte_carlo/include/criteria.h"

namespace feasst {

/**
  Metropolis acceptance criteria.
 */
class CriteriaMetropolis : public Criteria {
 public:
  CriteriaMetropolis(const argtype &args = argtype()) : Criteria(args) {}

  bool is_accepted(const Acceptance& acceptance,
    const System * system,
    const double uniform_random) override;

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<CriteriaMetropolis>(istr); }

  void serialize(std::ostream& ostr) const override;
  CriteriaMetropolis(std::istream& istr);
  ~CriteriaMetropolis() {}

 private:
  const std::string class_name_ = "CriteriaMetropolis";
};

inline std::shared_ptr<CriteriaMetropolis> MakeCriteriaMetropolis(const argtype &args = argtype()) {
  return std::make_shared<CriteriaMetropolis>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_CRITERIA_METROPOLIS_H_
