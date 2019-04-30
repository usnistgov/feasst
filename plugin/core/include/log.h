
#ifndef FEASST_CORE_LOG_H_
#define FEASST_CORE_LOG_H_

#include "core/include/analyze.h"

namespace feasst {

class Log : public AnalyzeWriteOnly {
 public:
  Log(const argtype &args = argtype()) : AnalyzeWriteOnly(args) {
    set_append();
  }
  void initialize(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    std::stringstream ss;
    ss << "#" << criteria->status_header() << " " << trial_factory.status_header()
       << std::endl;
    printer(ss.str());
  }

  std::string write(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    // ensure the following order matches the header from initialization.
    std::stringstream ss;
    ss << criteria->status() << " " << trial_factory.status() << std::endl;
    return ss.str();
  }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(1, ostr);
  }

  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    feasst_deserialize_version(istr);
    auto analyze = std::make_shared<Log>();
    return analyze;
  }

 private:
  const std::string class_name_ = "Log";
};

inline std::shared_ptr<Log> MakeLog(const argtype &args = argtype()) {
  return std::make_shared<Log>(args);
}

}  // namespace feasst

#endif  // FEASST_CORE_LOG_H_
