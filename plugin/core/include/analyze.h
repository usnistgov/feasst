
#ifndef FEASST_CORE_ANALYZE_H_
#define FEASST_CORE_ANALYZE_H_

#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include "core/include/trial_factory.h"
#include "core/include/trial_transfer.h"
#include "core/include/bond_visitor.h"
#include "core/include/file_xyz.h"
#include "core/include/stepper.h"
#include "core/include/arguments.h"

namespace feasst {

class Analyze : public Stepper {
 public:
  Analyze(const argtype &args = argtype()) : Stepper(args) {}

  virtual void initialize(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) {
    // do nothing by default
  }

  virtual void trial(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory);

  virtual void update(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory);

  virtual std::string write(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory);

  virtual ~Analyze() {}
};

class AnalyzeFactory : public Analyze {
 public:
  AnalyzeFactory() : Analyze() {}

  void initialize(const std::shared_ptr<Criteria> criteria,
    const System& system,
    const TrialFactory& trial_factory) override;

  void add(std::shared_ptr<Analyze> analyze) {
    analyzers_.push_back(analyze); }

  const std::vector<std::shared_ptr<Analyze> >& analyzers() const {
    return analyzers_; }

  void trial(const std::shared_ptr<Criteria> criteria,
    const System& system,
    const TrialFactory& trial_factory) override;

  private:
    std::vector<std::shared_ptr<Analyze> > analyzers_;
};

class AnalyzeWriteOnly : public Analyze {
 public:
  AnalyzeWriteOnly(
    /**
      steps_per : write every this many steps
     */
    const argtype &args = argtype()) : Analyze(args) {
    // disable update
    Analyze::set_steps_per_update(-1);

    // parse
    if (!args_.key("steps_per").empty()) {
      set_steps_per(args_.integer());
    }
  }

  void set_steps_per_update(const int steps) override {
    ERROR("This analyze is write only."); }

  void set_steps_per(const int steps) { set_steps_per_write(steps); }
};

class AnalyzeUpdateOnly : public Analyze {
 public:
  AnalyzeUpdateOnly(
    /**
      steps_per : update every this many steps
     */
    const argtype &args = argtype()) : Analyze(args) {
    // disable update
    Analyze::set_steps_per_write(-1);

    // parse
    if (!args_.key("steps_per").empty()) {
      set_steps_per(args_.integer());
    }
  }

  void set_steps_per_write(const int steps) override {
    ERROR("This analyze is update only."); }

  void set_steps_per(const int steps) { set_steps_per_update(steps); }
};

class Log : public AnalyzeWriteOnly {
 public:
  Log(const argtype &args = argtype()) : AnalyzeWriteOnly(args) {}
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
};

inline std::shared_ptr<Log> MakeLog(const argtype &args = argtype()) {
  return std::make_shared<Log>(args);
}

class Movie : public AnalyzeWriteOnly {
 public:
  Movie(const argtype &args = argtype()) : AnalyzeWriteOnly(args) {}
  void initialize(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    ASSERT(!file_name().empty(), "file name required. Did you forget to " <<
      "Analyze::set_file_name()?");
    xyz_.write(file_name(), system.configuration());
    xyz_.set_append(1);

    // write vmd
    std::stringstream ss;
    ss << file_name() << ".vmd";
    vmd_.write(ss.str(), system.configuration(), file_name());
  }

  std::string write(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    // ensure the following order matches the header from initialization.
    xyz_.write(file_name(), system.configuration());
    return std::string("");
  }

 private:
  FileXYZ xyz_;
  FileVMD vmd_;
};

inline std::shared_ptr<Movie> MakeMovie(const argtype &args = argtype()) {
  return std::make_shared<Movie>(args);
}

}  // namespace feasst

#endif  // FEASST_CORE_ANALYZE_H_
