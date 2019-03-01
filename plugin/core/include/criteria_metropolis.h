
#ifndef FEASST_CORE_CRITERIA_METROPOLIS_H_
#define FEASST_CORE_CRITERIA_METROPOLIS_H_

#include "core/include/criteria.h"
#include "core/include/random.h"

namespace feasst {

/**
  Metropolis acceptance criteria.
 */
class CriteriaMetropolis : public Criteria {
 public:
  CriteriaMetropolis(const argtype &args = argtype()) : Criteria(args) {}
  bool is_accepted(const AcceptanceCriteria accept_criteria) override;

 private:
  Random random_;
};

inline std::shared_ptr<CriteriaMetropolis> MakeCriteriaMetropolis(const argtype &args = argtype()) {
  return std::make_shared<CriteriaMetropolis>(args);
}

}  // namespace feasst

#endif  // FEASST_CORE_CRITERIA_METROPOLIS_H_
