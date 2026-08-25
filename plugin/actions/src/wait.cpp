#include <chrono>
#include <thread>
#include "utils/include/file.h"  // file_exists
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "actions/include/wait.h"

namespace feasst {

Wait::Wait(argtype * args) {
  class_name_ = "Wait";
  until_file_exists_ = str("until_file_exists", args, "");
  seconds_per_check_ = dble("seconds_per_check", args, 5);
  std::vector<std::string> sdat = split(until_file_exists_, '[');
  if (static_cast<int>(sdat.size()) == 2) {
    std::vector<std::string> sdat2 = split(sdat[1], ']');
    if (static_cast<int>(sdat2.size()) == 2) {
      prefix_ = sdat[0];
      suffix_ = sdat2[1];
      number_ = str_to_int(sdat2[0]);
      ASSERT(number_ > 0, "number: " << number_);
    }
    DEBUG("sdat2 " << feasst_str(sdat2));
  }
  DEBUG("sdat " << feasst_str(sdat));
}
Wait::Wait(argtype args) : Wait(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(Wait,);

Wait::Wait(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 6345 && version <= 6345, "mismatch version: " << version);
  feasst_deserialize(&until_file_exists_, istr);
  feasst_deserialize(&seconds_per_check_, istr);
  feasst_deserialize(&prefix_, istr);
  feasst_deserialize(&suffix_, istr);
  feasst_deserialize(&number_, istr);
}

void Wait::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(6345, ostr);
  feasst_serialize(until_file_exists_, ostr);
  feasst_serialize(seconds_per_check_, ostr);
  feasst_serialize(prefix_, ostr);
  feasst_serialize(suffix_, ostr);
  feasst_serialize(number_, ostr);
}

void Wait::run(MonteCarlo * mc) {
  std::chrono::duration<double> duration(seconds_per_check_);
  if (prefix_.empty()) {
    if (!file_exists(until_file_exists_)) {
      while (!file_exists(until_file_exists_)) {
        std::this_thread::sleep_for(duration);
      }
      std::this_thread::sleep_for(duration);
    }
  } else {
    bool any_false = false;
    for (int i = 0; i < number_; ++i) {
      const std::string filename = prefix_ + str(i) + suffix_;
      while (!file_exists(filename)) {
        std::this_thread::sleep_for(duration);
        any_false = true;
      }
      DEBUG("found " << filename);
    }
    if (any_false) {
      std::this_thread::sleep_for(duration);
    }
  }
}

}  // namespace feasst
