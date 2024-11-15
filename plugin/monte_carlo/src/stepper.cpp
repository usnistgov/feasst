#include <string>
#include <iostream>
#include <fstream>
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "math/include/accumulator.h"
#include "system/include/system.h"
#include "monte_carlo/include/stepper.h"
#include "monte_carlo/include/criteria.h"

namespace feasst {

Stepper::Stepper(argtype args) : Stepper(&args) {
  feasst_check_all_used(args);
}
Stepper::Stepper(argtype * args) {
  set_trials_per_write(integer("trials_per_write", args, 1));
  set_trials_per_update(integer("trials_per_update", args, 1));
  if (used("output_file", *args)) {
    set_output_file(str("output_file", args));
  }
  if (used("file_name", *args)) {
    WARN("Stepper argument file_name was renamed to output_file.");
    set_output_file(str("file_name", args));
  }

  if (boolean("append", args, false)) {
    set_append();
  } else {
    set_no_append();
  }

  if (boolean("clear_file", args, false)) {
    ASSERT(!output_file_.empty(), "output_file is a required argument with clear.");
    std::ofstream file;
    file.open(output_file_, std::ofstream::out);
    file.close();
  }
  rewrite_header_ = boolean("rewrite_header", args, true);
  stop_after_phase_ = integer("stop_after_phase", args, -1);
  start_after_phase_ = integer("start_after_phase", args, -1);
  stop_after_iteration_ = integer("stop_after_iteration", args, -1);
  start_after_iteration_ = integer("start_after_iteration", args, -1);
  output_file_append_phase_ =
    boolean("output_file_append_phase", args, false);
  if (used("file_name_append_phase", *args)) {
    WARN("Stepper::file_name_append_phase was renamed to output_file_append_phase.");
    output_file_append_phase_ = boolean("file_name_append_phase", args);
  }
  set_multistate(boolean("multistate", args, false));
  is_multistate_aggregate_ =
    boolean("multistate_aggregate", args, true);
  if (is_multistate() && is_multistate_aggregate_) {
    rewrite_header_ = false;
  }
  accumulator_ = std::make_shared<Accumulator>(args);
//  if (used("block_size", *args)) {
//    accumulator_.set_block_size(integer("block_size", args));
//  }
//  if (used("num_moments", *args)) {
//    accumulator_.set_moments(integer("num_moments", args));
//  }
  configuration_index_ = integer("configuration_index", args, 0);
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

std::string Stepper::output_file(const Criteria& criteria) const {
  if (output_file_.empty() || !output_file_append_phase_) {
    return output_file_;
  }
  return output_file_ + "_phase" + str(criteria.phase());
}

void Stepper::printer(const std::string output, const std::string& output_file) {
  DEBUG("output_file " << output_file);
  if (output_file.empty() && !is_multistate_aggregate()) {
    std::cout << output;
  } else {
    std::ofstream file;
    if (append_) {
      file.open(output_file, std::ofstream::out | std::ofstream::app);
    } else {
      file.open(output_file, std::ofstream::out);
    }
    file << output;
    file.close();
  }
}

void Stepper::set_state(const int state) {
  state_ = state;
  if (!output_file_.empty()) {
    std::stringstream ss;
    ss << output_file_ << "_state" << state;
    output_file_ = ss.str();
  }
}

void Stepper::serialize(std::ostream& ostr) const {
  ostr << class_name() << " ";
  feasst_serialize_version(497, ostr);
  feasst_serialize(trials_since_update_, ostr);
  feasst_serialize(trials_since_write_, ostr);
  feasst_serialize(trials_per_update_, ostr);
  feasst_serialize(trials_per_write_, ostr);
  feasst_serialize(output_file_, ostr);
  feasst_serialize(append_, ostr);
  feasst_serialize(stop_after_phase_, ostr);
  feasst_serialize(start_after_phase_, ostr);
  feasst_serialize(stop_after_iteration_, ostr);
  feasst_serialize(start_after_iteration_, ostr);
  feasst_serialize(output_file_append_phase_, ostr);
  feasst_serialize(is_multistate_, ostr);
  feasst_serialize(is_multistate_aggregate_, ostr);
  feasst_serialize(state_, ostr);
  feasst_serialize(configuration_index_, ostr);
  feasst_serialize(rewrite_header_, ostr);
  feasst_serialize(accumulator_, ostr);
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
  feasst_deserialize(&output_file_, istr);
  feasst_deserialize(&append_, istr);
  feasst_deserialize(&stop_after_phase_, istr);
  feasst_deserialize(&start_after_phase_, istr);
  feasst_deserialize(&stop_after_iteration_, istr);
  feasst_deserialize(&start_after_iteration_, istr);
  feasst_deserialize(&output_file_append_phase_, istr);
  feasst_deserialize(&is_multistate_, istr);
  feasst_deserialize(&is_multistate_aggregate_, istr);
  feasst_deserialize(&state_, istr);
  feasst_deserialize(&configuration_index_, istr);
  feasst_deserialize(&rewrite_header_, istr);
//  feasst_deserialize(accumulator_, istr);
// HWH for unknown reasons, this function template does not work.
  {
    int existing;
    istr >> existing;
    if (existing != 0) {
      accumulator_ = std::make_shared<Accumulator>(istr);
    }
  }
  feasst_deserialize_endcap("Stepper", istr);
}

const Configuration& Stepper::configuration(const System& system) const {
  return system.configuration(configuration_index_);
}

const Accumulator& Stepper::accumulator() const { return *accumulator_; }

Accumulator * Stepper::get_accumulator() { return accumulator_.get(); }

}  // namespace feasst
