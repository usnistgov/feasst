#include <iostream>
#include <iomanip>
#include <fstream>
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "utils/include/progress_report.h"

namespace feasst {

void ProgressReport::reset() {
  last_percent_ = -1.;
  current_ = 0.;
}

ProgressReport::ProgressReport(argtype args) {
  reset();
  num_ = integer("num", &args, 1);
  percent_per_write_ = dble("percent_per_write", &args, 0.1);
  if (used("file_name", args)) file_name_ = str("file_name", &args);
  ASSERT(percent_per_write_ > 0.0001,
    "percent_per_write: " << percent_per_write_ << " must be > 0.0001");
  FEASST_CHECK_ALL_USED(args);
}

void ProgressReport::serialize(std::ostream& ostr) const {
  feasst_serialize_version(5756, ostr);
  feasst_serialize(num_, ostr);
  feasst_serialize(current_, ostr);
  feasst_serialize(percent_per_write_, ostr);
  feasst_serialize(file_name_, ostr);
  feasst_serialize(last_percent_, ostr);
}

ProgressReport::ProgressReport(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 5756, "version mismatch: " << version);
  feasst_deserialize(&num_, istr);
  feasst_deserialize(&current_, istr);
  feasst_deserialize(&percent_per_write_, istr);
  feasst_deserialize(&file_name_, istr);
  feasst_deserialize(&last_percent_, istr);
}

double ProgressReport::percent() const {
  return static_cast<double>(current_)/static_cast<double>(num_);
}

void ProgressReport::check() {
  if (current_ == 0) {
    starting_hours_ = cpu_hours();
  }
  ++current_;
  if (current_ == num_ || percent() >= last_percent_ + percent_per_write_) {
    write();
    last_percent_ = percent();
  }
}

void ProgressReport::write() {
  std::stringstream ss;
  const double elapsed_hours = cpu_hours() - starting_hours_;
  if (current_ == 1) {
    ss << "#percent,hours_elapsed,hours_remain";
  } else if (current_ == num_) {
    ss << "1 " << elapsed_hours << " 0";
  } else {
    DEBUG("elapsed " << elapsed_hours);
    const double percent_per_hours = percent()/elapsed_hours;
    DEBUG("percent_per_hours " << percent_per_hours);
    const double remaining_hours = (1. - percent())/percent_per_hours;
    ss << std::setprecision(3) <<  percent() << " "
       << elapsed_hours << " " << remaining_hours;
  }
  DEBUG("filename? " << file_name_);
  if (file_name_.empty()) {
    std::cout << ss.str() << std::endl;
  } else {
    std::ofstream file;
    file.open(file_name_, std::ofstream::out | std::ofstream::app);
    file << ss.str() << std::endl;
    file.close();
  }
}

}  // namespace feasst
