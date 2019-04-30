
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
  CriteriaMayer(const argtype& args = argtype()) : Criteria(args) {}

  bool is_accepted(const AcceptanceCriteria accept_criteria) override;

  double second_virial() const {
    return (2./3.)*PI*mayer_.average()/mayer_ref_.average();
  }

  bool verbose = false;

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<CriteriaMayer>(istr); }
  void serialize(std::ostream& ostr) const override;
  CriteriaMayer(std::istream& istr);
  ~CriteriaMayer() {}

 private:
  const std::string class_name_ = "CriteriaMayer";
  double f12old_ = -1.;
  double f12ref_ = -1.;
  Accumulator mayer_;
  Accumulator mayer_ref_;
  int reference_index_ = 0;

  // not currently serialized
  Random random_;
};

inline std::shared_ptr<CriteriaMayer> MakeCriteriaMayer(const argtype &args = argtype()) {
  return std::make_shared<CriteriaMayer>(args);
}

}  // namespace feasst

#endif  // FEASST_CORE_CRITERIA_MAYER_H_
