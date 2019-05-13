
#ifndef FEASST_STEPPERS_CRITERIA_WRITER_H_
#define FEASST_STEPPERS_CRITERIA_WRITER_H_

#include "monte_carlo/include/analyze.h"

namespace feasst {

class CriteriaWriter : public AnalyzeWriteOnly {
 public:
  CriteriaWriter(const argtype &args = argtype());
  std::string write(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    // ensure the following order matches the header from initialization.
    std::stringstream ss;
    ss << criteria->write() << std::endl;
    return ss.str();
  }

  void serialize(std::ostream& ostr) const override {
    Stepper::serialize(ostr);
    feasst_serialize_version(1, ostr);
  }

  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<CriteriaWriter>(istr); }

  CriteriaWriter(std::istream& istr) : AnalyzeWriteOnly(istr) {
    feasst_deserialize_version(istr); }

  std::string class_name() const override { return std::string("CriteriaWriter"); }
};

inline std::shared_ptr<CriteriaWriter> MakeCriteriaWriter(const argtype &args = argtype()) {
  return std::make_shared<CriteriaWriter>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_CRITERIA_WRITER_H_
