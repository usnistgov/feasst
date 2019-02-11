
#ifndef FEASST_CORE_ANALYZE_H_
#define FEASST_CORE_ANALYZE_H_

#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include "core/include/trial_factory.h"
#include "core/include/trial_transfer.h"

namespace feasst {

class Analyze {
 public:
  Analyze() { set_steps_per_update(); }

  virtual void initialize(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) {
    // do nothing by default
  }

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

  /// Write to standard output if file name is not set. Otherwise, to file.
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

  virtual void trial(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) {
    if (is_time(steps_per_update_, &steps_since_update_)) {
      update(criteria, system, trial_factory);
    }
    if (is_time(steps_per_write_, &steps_since_write_)) {
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

  /// Set the number of trial steps per analysis update.
  /// Disabled if steps is not positive.
  virtual void set_steps_per_update(const int steps = 1) { steps_per_update_ = steps; }

  /// Set the number of trial steps per writing of analysis to file or screen.
  /// Disabled if steps is not positive.
  void set_steps_per_write(const int steps = 1) { steps_per_write_ = steps; }

  /// Set the name of the file to write. If empty, write to screen.
  void set_file_name(const std::string file_name) { file_name_ = file_name; }

  virtual ~Analyze() {}

 private:
  int steps_per_update_;
  int steps_per_write_;
  int steps_since_update_ = 0;
  int steps_since_write_ = 0;
  std::string file_name_;
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
};

class Log : public AnalyzeWriteOnly {
 public:
  void initialize(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    std::stringstream ss;
    ss << "#" << system.status_header() << trial_factory.status_header()
       << std::endl;
    printer(ss.str());
  }

  std::string write(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    // ensure the following order matches the header from initialization.
    std::stringstream ss;
    ss << system.status() << trial_factory.status() << std::endl;
    return ss.str();
  }
};

}  // namespace feasst

#endif  // FEASST_CORE_ANALYZE_H_
