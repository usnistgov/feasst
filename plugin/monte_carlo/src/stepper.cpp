#include <string>
#include <fstream>
#include <sstream>
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "math/include/accumulator.h"
#include "monte_carlo/include/stepper.h"
#include "monte_carlo/include/criteria.h"

namespace feasst {

Stepper::Stepper(argtype args) : Stepper(&args) {
  FEASST_CHECK_ALL_USED(args);
}
Stepper::Stepper(argtype * args) {
  set_trials_per_write(integer("trials_per_write", args, 1));
  set_trials_per_update(integer("trials_per_update", args, 1));
  if (used("file_name", *args)) set_file_name(str("file_name", args));

  if (boolean("append", args, false)) {
    set_append();
  } else {
    set_no_append();
  }

  if (boolean("clear_file", args, false)) {
    ASSERT(!file_name_.empty(), "file_name is a required argument with clear.");
    std::ofstream file;
    file.open(file_name_, std::ofstream::out);
    file.close();
  }
  rewrite_header_ = boolean("rewrite_header", args, true);
  stop_after_phase_ = integer("stop_after_phase", args, -1);
  start_after_phase_ = integer("start_after_phase", args, -1);
  file_name_append_phase_ =
    boolean("file_name_append_phase", args, false);
  set_multistate(boolean("multistate", args, false));
  is_multistate_aggregate_ =
    boolean("multistate_aggregate", args, true);
  if (is_multistate() && is_multistate_aggregate_) {
    rewrite_header_ = false;
  }
  accumulator_ = Accumulator(args);
//  if (used("block_size", *args)) {
//    accumulator_.set_block_size(integer("block_size", args));
//  }
//  if (used("num_moments", *args)) {
//    accumulator_.set_moments(integer("num_moments", args));
//  }
  configuration_ = integer("configuration", args, 0);
}

bool Stepper::is_time(const int trials_per, int * trials_since) {
  if (trials_per > 0) {
    ++(*trials_since);
    if (*trials_since >= trials_per) {
      *trials_since = 0;
      return true;
    }
  }
  return false;
}

std::string Stepper::file_name(const Criteria& criteria) const {
  if (file_name_.empty() || !file_name_append_phase_) {
    return file_name_;
  }
  return file_name_ + "_phase" + str(criteria.phase());
}

void Stepper::printer(const std::string output, const std::string& file_name) {
  DEBUG("filename? " << file_name);
  if (file_name.empty() && !is_multistate_aggregate()) {
    std::cout << output;
  } else {
    std::ofstream file;
    if (append_) {
      file.open(file_name, std::ofstream::out | std::ofstream::app);
    } else {
      file.open(file_name, std::ofstream::out);
    }
    file << output;
    file.close();
  }
}

void Stepper::set_state(const int state) {
  state_ = state;
  if (!file_name_.empty()) {
    std::stringstream ss;
    ss << file_name_ << "_state" << state;
    file_name_ = ss.str();
  }
}

void Stepper::serialize(std::ostream& ostr) const {
  ostr << class_name() << " ";
  feasst_serialize_version(497, ostr);
  feasst_serialize(trials_since_update_, ostr);
  feasst_serialize(trials_since_write_, ostr);
  feasst_serialize(trials_per_update_, ostr);
  feasst_serialize(trials_per_write_, ostr);
  feasst_serialize(file_name_, ostr);
  feasst_serialize(append_, ostr);
  feasst_serialize(stop_after_phase_, ostr);
  feasst_serialize(start_after_phase_, ostr);
  feasst_serialize(file_name_append_phase_, ostr);
  feasst_serialize(is_multistate_, ostr);
  feasst_serialize(is_multistate_aggregate_, ostr);
  feasst_serialize(state_, ostr);
  feasst_serialize(configuration_, ostr);
  feasst_serialize(rewrite_header_, ostr);
  feasst_serialize_fstobj(accumulator_, ostr);
  feasst_serialize_endcap("Stepper", ostr);
}

Stepper::Stepper(std::istream& istr) {
  std::string name;
  istr >> name;
  const int version = feasst_deserialize_version(istr);
  ASSERT(497 == version, version);
  feasst_deserialize(&trials_since_update_, istr);
  feasst_deserialize(&trials_since_write_, istr);
  feasst_deserialize(&trials_per_update_, istr);
  feasst_deserialize(&trials_per_write_, istr);
  feasst_deserialize(&file_name_, istr);
  feasst_deserialize(&append_, istr);
  feasst_deserialize(&stop_after_phase_, istr);
  feasst_deserialize(&start_after_phase_, istr);
  feasst_deserialize(&file_name_append_phase_, istr);
  feasst_deserialize(&is_multistate_, istr);
  feasst_deserialize(&is_multistate_aggregate_, istr);
  feasst_deserialize(&state_, istr);
  feasst_deserialize(&configuration_, istr);
  feasst_deserialize(&rewrite_header_, istr);
  feasst_deserialize_fstobj(&accumulator_, istr);
  feasst_deserialize_endcap("Stepper", istr);
}

}  // namespace feasst
