
#ifndef FEASST_CORE_MODIFY_H_
#define FEASST_CORE_MODIFY_H_

#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include <time.h>
#include "core/include/stepper.h"
#include "core/include/criteria.h"
#include "core/include/trial_factory.h"

namespace feasst {

class Modify : public Stepper {
 public:
  Modify(const argtype &args = argtype()) : Stepper(args) {}

  virtual void initialize(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) {
    // do nothing by default
  }

  virtual void trial(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) {
    if (is_time(steps_per_update(), &steps_since_update_)) {
      update(criteria, system, trial_factory);
    }
    if (is_time(steps_per_write(), &steps_since_write_)) {
      printer(write(criteria, system, trial_factory));
    }
  }

  virtual void update(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) {
    ERROR("not implemented");
  }

  virtual std::string write(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) {
    ERROR("not implemented");
    return std::string("");
  }

  std::string class_name() const override { return std::string("Modify"); }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Modify> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Modify> >& deserialize_map();
  std::shared_ptr<Modify> deserialize(std::istream& istr);
  Modify(std::istream& istr) : Stepper(istr) {}
  virtual ~Modify() {}
};

class ModifyUpdateOnly : public Modify {
 public:
  ModifyUpdateOnly(
    /**
      steps_per : update every this many steps
     */
    const argtype &args = argtype()) : Modify(args) {
    // disable write
    Modify::set_steps_per_write(-1);

    // parse
    if (!args_.key("steps_per").empty()) {
      set_steps_per(args_.integer());
    }
  }

  void set_steps_per_write(const int steps) override {
    ERROR("This modify is update only."); }

  void set_steps_per(const int steps) { set_steps_per_update(steps); }

  ModifyUpdateOnly(std::istream& istr) : Modify(istr) {}
};

}  // namespace feasst

#endif  // FEASST_CORE_MODIFY_H_
