
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
  bool is_accepted(const AcceptanceCriteria accept_criteria) override;

 private:
  Random random_;
};

}  // namespace feasst

#endif  // FEASST_CORE_CRITERIA_METROPOLIS_H_
