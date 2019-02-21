
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
#include "core/include/debug.h"

namespace feasst {

class Stepper {
 public:
  Stepper() { set_steps_per_update(); }

  /// Check if it is time to update or write. Return true if so.
  bool is_time(const int steps_per, int * steps_since) {
    if (steps_per > 0) {
      ++(*steps_since);
      if (*steps_since == steps_per) {
        *steps_since = 0;
        return true;
      } else {
        ASSERT(*steps_since < steps_per,
          "skipped an analysis step?");
      }
    }
    return false;
  }

  /// Write to standard output if file name is not set. Otherwise, output file.
  void printer(const std::string output) {
    if (file_name_.empty()) {
      std::cout << output;
    } else {
      std::ofstream file;
      file.open(file_name_, std::ofstream::out | std::ofstream::app);
      file << output;
      file.close();
    }
  }

  /// Set the number of trial steps per analysis update.
  /// Disabled if steps is not positive.
  virtual void set_steps_per_update(const int steps = 1) { steps_per_update_ = steps; }

  /// Set the number of trial steps per writing of analysis to file or screen.
  /// Disabled if steps is not positive.
  virtual void set_steps_per_write(const int steps = 1) { steps_per_write_ = steps; }

  /// Set the name of the file to write. If empty, write to screen.
  void set_file_name(const std::string file_name) { file_name_ = file_name; }
  std::string file_name() const { return file_name_; }

  int steps_per_update() const { return steps_per_update_; }
  int steps_per_write() const { return steps_per_write_; }
  int * get_steps_since_update() { return &steps_since_update_; }
  int * get_steps_since_write() { return &steps_since_write_; }

  virtual ~Stepper() {}

 private:
  int steps_per_update_;
  int steps_per_write_;
  int steps_since_update_ = 0;
  int steps_since_write_ = 0;
  std::string file_name_;
};

class Analyze : public Stepper {
 public:
  virtual void initialize(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) {
    // do nothing by default
  }

  virtual void trial(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) {
    if (is_time(steps_per_update(), get_steps_since_update())) {
      update(criteria, system, trial_factory);
    }
    if (is_time(steps_per_write(), get_steps_since_write())) {
      printer(write(criteria, system, trial_factory));
    }
  }

  virtual void update(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) {
    ERROR("not implemented");
  }

  virtual std::string write(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) {
    ERROR("not implemented");
    return std::string("");
  }

  virtual ~Analyze() {}
};

class AnalyzeFactory : public Analyze {
 public:
  AnalyzeFactory() : Analyze() {}

  void add(std::shared_ptr<Analyze> analyze) {
    analyzers_.push_back(analyze);
  }

  const std::vector<std::shared_ptr<Analyze> >& analyzers() const {
    return analyzers_; }

  void trial(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    for (const std::shared_ptr<Analyze> analyze : analyzers_) {
      analyze->trial(criteria, system, trial_factory);
    }
  }

  private:
    std::vector<std::shared_ptr<Analyze> > analyzers_;
};

class AnalyzeWriteOnly : public Analyze {
 public:
  AnalyzeWriteOnly() : Analyze() {
    // disable update
    Analyze::set_steps_per_update(-1); }

  void set_steps_per_update(const int steps) override {
    ERROR("This analyze is write only."); }

  void set_steps_per(const int steps) { set_steps_per_write(steps); }
};

class AnalyzeUpdateOnly : public Analyze {
 public:
  AnalyzeUpdateOnly() : Analyze() {
    // disable update
    Analyze::set_steps_per_write(-1); }

  void set_steps_per_write(const int steps) override {
    ERROR("This analyze is update only."); }

  void set_steps_per(const int steps) { set_steps_per_update(steps); }
};

class Log : public AnalyzeWriteOnly {
 public:
  void initialize(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    std::stringstream ss;
    ss << "#" << criteria->status_header() << " " << system.status_header() << " " << trial_factory.status_header()
       << std::endl;
    printer(ss.str());
  }

  std::string write(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    // ensure the following order matches the header from initialization.
    std::stringstream ss;
    ss << criteria->status() << " " << system.status() << " " << trial_factory.status() << std::endl;
    return ss.str();
  }
};

class Movie : public AnalyzeWriteOnly {
 public:
  void initialize(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    ASSERT(!file_name().empty(), "file name required");
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

class RigidBondChecker : public AnalyzeUpdateOnly {
 public:
  void update(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    visitor_.compute(bond_, system.configuration());
    ASSERT(std::abs(visitor_.energy()) < NEAR_ZERO, "bond check failure");
    visitor_.compute(angle_, system.configuration());
    ASSERT(std::abs(visitor_.energy()) < NEAR_ZERO, "angle check failure");
  }

 private:
  BondVisitor visitor_;
  BondSquareWell bond_;
  AngleSquareWell angle_;
};

}  // namespace feasst

#endif  // FEASST_CORE_ANALYZE_H_
