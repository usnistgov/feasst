
#ifndef FEASST_CORE_STEPPER_H_
#define FEASST_CORE_STEPPER_H_

#include <string>
#include "core/include/arguments.h"

namespace feasst {

/**
  Perform an action (update or write) every so many steps.
  This action could be read-only (see analyze) or not (see modify).
  Write to screen if file name is not provided.
 */
class Stepper {
 public:
  Stepper(
    /**
      steps_per_write : write every this many steps
      steps_per_update : update every this many steps
      file_name : file name to save output
      append : append file output if set to 1. Do not append if 0 (default).
     */
    const argtype &args = argtype());

  /// Check if it is time to update or write. Return true if so.
  bool is_time(const int steps_per, int * steps_since);

  /// Write to standard output if file name is not set. Otherwise, output file.
  void printer(const std::string output);

  /// Set the number of trial steps per analysis update.
  /// Disabled if steps is not positive.
  virtual void set_steps_per_update(const int steps = 1) {
    steps_per_update_ = steps; }

  /// Set the number of trial steps per writing of analysis to file or screen.
  /// Disabled if steps is not positive.
  virtual void set_steps_per_write(const int steps) {
    steps_per_write_ = steps; }

  /// Set the name of the file to write. If empty, write to screen.
  void set_file_name(const std::string file_name) { file_name_ = file_name; }
  std::string file_name() const { return file_name_; }

  /// Return the number of steps per update
  int steps_per_update() const { return steps_per_update_; }

  /// Return the number of steps per write.
  int steps_per_write() const { return steps_per_write_; }

  /// Set file output to append.
  void set_append() { append_ = true; }

  /// Set file output to not append.
  void set_no_append() { append_ = false; }

  void serialize(std::ostream& ostr) const;
  Stepper(std::istream& istr);
  virtual ~Stepper() {}

 protected:
  Arguments args_;
  int steps_since_update_ = 0;
  int steps_since_write_ = 0;

 private:
  int steps_per_update_;
  int steps_per_write_;
  std::string file_name_;
  bool append_;
};

}  // namespace feasst

#endif  // FEASST_CORE_STEPPER_H_
