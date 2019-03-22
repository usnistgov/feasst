#include <string>
#include <fstream>
#include "core/include/stepper.h"
#include "core/include/debug.h"

namespace feasst {

Stepper::Stepper(const argtype &args) {
  // defaults
  set_steps_per_update();
  set_no_append();

  // parse
  args_.init(args);
  if (args_.key("steps_per_write").used()) {
    set_steps_per_write(args_.integer());
  }
  if (args_.key("steps_per_update").used()) {
    set_steps_per_update(args_.integer());
  }
  if (args_.key("file_name").used()) {
    set_file_name(args_.str());
  }
  if (args_.key("append").used()) {
    const int flag = args_.integer();
    if (flag == 0) {
      set_no_append();
    } else if (flag == 1) {
      set_append();
    } else {
      ERROR("unrecognized append argument: " << flag);
    }
  }
}

bool Stepper::is_time(const int steps_per, int * steps_since) {
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

void Stepper::printer(const std::string output) {
  if (file_name_.empty()) {
    std::cout << output;
  } else {
    std::ofstream file;
    if (append_) {
      file.open(file_name_, std::ofstream::out | std::ofstream::app);
    } else {
      file.open(file_name_, std::ofstream::out);
    }
    file << output;
    file.close();
  }
}

}  // namespace feasst
