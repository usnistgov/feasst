
#ifndef FEASST_CORE_CRITERIA_MAYER_H_
#define FEASST_CORE_CRITERIA_MAYER_H_

#include "core/include/criteria.h"
#include "core/include/random.h"
#include "core/include/accumulator.h"
#include "core/include/model_hard_sphere.h"

namespace feasst {

/**
  Mayer-sampling Monte Carlo acceptance criteria.
  Currently assumes hard sphere reference.
  HWH add reference.
 */
class CriteriaMayer : public Criteria {
 public:
  bool is_accepted(const AcceptanceCriteria accept_criteria) override;

  double second_virial() const {
    return (2./3.)*PI*mayer_.average()/mayer_ref_.average();
  }

  bool verbose = false;

  ~CriteriaMayer() {}
 private:
  Random random_;
  double f12old_ = -1.;
  double f12ref_ = -1.;
  Accumulator mayer_;
  Accumulator mayer_ref_;
  int reference_index_ = 0;
};

}  // namespace feasst

#endif  // FEASST_CORE_CRITERIA_MAYER_H_
