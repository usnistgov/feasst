
#ifndef FEASST_MAYER_CRITERIA_MAYER_H_
#define FEASST_MAYER_CRITERIA_MAYER_H_

#include "monte_carlo/include/criteria.h"
#include "math/include/accumulator.h"
#include "system/include/model_hard_sphere.h"

namespace feasst {

/**
  Mayer-sampling Monte Carlo acceptance criteria.
  Currently assumes hard sphere reference.
  HWH add reference.
 */
class CriteriaMayer : public Criteria {
 public:
  CriteriaMayer(const argtype& args = argtype()) : Criteria(args) {}

  bool is_accepted(const Acceptance& acceptance,
    const System * system,
    const double uniform_random) override;

  double second_virial() const {
    return (2./3.)*PI*mayer_.average()/mayer_ref_.average();
  }

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
};

inline std::shared_ptr<CriteriaMayer> MakeCriteriaMayer(const argtype &args = argtype()) {
  return std::make_shared<CriteriaMayer>(args);
}

}  // namespace feasst

#endif  // FEASST_MAYER_CRITERIA_MAYER_H_
