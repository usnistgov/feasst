
#include "core/include/debug.h"
#include "core/include/analyze.h"

namespace feasst {

std::map<std::string, std::shared_ptr<Analyze> >& Analyze::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Analyze> >* ans =
     new std::map<std::string, std::shared_ptr<Analyze> >();
  return *ans;
}

void Analyze::serialize(std::ostream& ostr) const { ERROR("not implemented"); }

std::shared_ptr<Analyze> Analyze::create(std::istream& istr) const {
  ERROR("not implemented");
}

std::shared_ptr<Analyze> Analyze::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr);
}

void Analyze::trial(const std::shared_ptr<Criteria> criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  if (is_time(steps_per_update(), &steps_since_update_)) {
    update(criteria, system, trial_factory);
  }
  if (is_time(steps_per_write(), &steps_since_write_)) {
    printer(write(criteria, system, trial_factory));
  }
}

void Analyze::update(const std::shared_ptr<Criteria> criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  ERROR("not implemented");
}

std::string Analyze::write(const std::shared_ptr<Criteria> criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  ERROR("not implemented");
  return std::string("");
}

class MapAnalyzeWriteOnly {
 public:
  MapAnalyzeWriteOnly() {
    AnalyzeWriteOnly().deserialize_map()["AnalyzeWriteOnly"] =
      std::make_shared<AnalyzeWriteOnly>();
  }
};

static MapAnalyzeWriteOnly mapper_ = MapAnalyzeWriteOnly();

std::shared_ptr<Analyze> AnalyzeWriteOnly::create(std::istream& istr) const {
  feasst_deserialize_version(istr);
  auto analyze = std::make_shared<AnalyzeWriteOnly>();
  return analyze;
}

void AnalyzeWriteOnly::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(1, ostr);
}

class MapAnalyzeUpdateOnly {
 public:
  MapAnalyzeUpdateOnly() {
    AnalyzeUpdateOnly().deserialize_map()["AnalyzeUpdateOnly"] =
      std::make_shared<AnalyzeUpdateOnly>();
  }
};

static MapAnalyzeUpdateOnly mapper_analyze_update_only = MapAnalyzeUpdateOnly();

std::shared_ptr<Analyze> AnalyzeUpdateOnly::create(std::istream& istr) const {
  feasst_deserialize_version(istr);
  auto model = std::make_shared<AnalyzeUpdateOnly>();
  return model;
}

void AnalyzeUpdateOnly::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(1, ostr);
}

}  // namespace feasst
