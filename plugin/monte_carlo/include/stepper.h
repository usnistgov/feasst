
#ifndef FEASST_MONTE_CARLO_STEPPER_H_
#define FEASST_MONTE_CARLO_STEPPER_H_

#include <string>
#include "utils/include/arguments.h"
#include "math/include/accumulator.h"

namespace feasst {

/**
  Perform an action (update or write) every so many steps.
  This action could be read-only (see analyze) or not (see modify).
  Write to screen if file name is not provided.
 */
class Stepper {
 public:
  /**
    args:
    - steps_per_write: Set the number of trial steps per write (default: 1).
      Disabled if negative value is provided.
    - steps_per_update: Set the number of trial steps per update (default: 1).
      Disabled if negative value is provided.
    - file_name: Set the file name to write output (default: empty).
      If empty, write to screen.
    - append: append file output if set to true.
      Do not append if false (default: "false").
    - clear_file: set true to clear contents of file_name, if exists.
      (default: false).
    - multistate:  set "true" to copy for each state (default: "false")
   */
  Stepper(const argtype &args = argtype());

  /// Return the number of steps per update
  int steps_per_update() const { return steps_per_update_; }

  /// Return the number of steps per write.
  int steps_per_write() const { return steps_per_write_; }

  /// Return the file name.
  std::string file_name() const { return file_name_; }

  /// Return true if appending.
  bool append() const { return append_; }

  /// Set the state. Append file name if not empty.
  void set_state(const int state = 0);

  /// Return if multistate.
  bool is_multistate() const { return is_multistate_; }

  /// Return the state.
  int state() const { return state_; }

  /// Return the accumulator.
  virtual const Accumulator& accumulator() const { FATAL("not implemented"); }

  virtual std::string class_name() const { return std::string("Stepper"); }

  void serialize(std::ostream& ostr) const;
  Stepper(std::istream& istr);
  virtual ~Stepper() {}

 protected:
  Arguments args_;
  int steps_since_update_ = 0;
  int steps_since_write_ = 0;

  /// Note that this should not be called after set_state, which appends name.
  void set_file_name(const std::string file_name) { file_name_ = file_name; }

  virtual void set_steps_per_update(const int steps = 1) {
    steps_per_update_ = steps; }

  virtual void set_steps_per_write(const int steps) {
    steps_per_write_ = steps; }

  /// Check if it is time to update or write. Return true if so.
  bool is_time(const int steps_per, int * steps_since);

  /// Write to standard output if file name is not set. Otherwise, output file.
  void printer(const std::string output);

  /// Set file output to append.
  void set_append() { append_ = true; }

  /// Set file output to not append.
  void set_no_append() { append_ = false; }

  /// Replicate the stepper individually for each state during initialization
  /// of the factory.
  void set_multistate(const bool multi) { is_multistate_ = multi; }

 private:
  int steps_per_update_;
  int steps_per_write_;
  std::string file_name_;
  bool append_;
  bool is_multistate_;
  int state_ = 0;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_STEPPER_H_
