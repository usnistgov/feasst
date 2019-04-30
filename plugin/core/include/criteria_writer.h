
#ifndef FEASST_CORE_CRITERIA_WRITER_H_
#define FEASST_CORE_CRITERIA_WRITER_H_

#include "core/include/analyze.h"

namespace feasst {

class CriteriaWriter : public AnalyzeWriteOnly {
 public:
  CriteriaWriter(const argtype &args = argtype()) : AnalyzeWriteOnly(args) {}
  std::string write(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    // ensure the following order matches the header from initialization.
    std::stringstream ss;
    ss << criteria->write() << std::endl;
    return ss.str();
  }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(1, ostr);
  }

  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    feasst_deserialize_version(istr);
    auto analyze = std::make_shared<CriteriaWriter>();
    return analyze;
  }

 private:
  const std::string class_name_ = "CriteriaWriter";
};

inline std::shared_ptr<CriteriaWriter> MakeCriteriaWriter(const argtype &args = argtype()) {
  return std::make_shared<CriteriaWriter>(args);
}

}  // namespace feasst

#endif  // FEASST_CORE_CRITERIA_WRITER_H_
